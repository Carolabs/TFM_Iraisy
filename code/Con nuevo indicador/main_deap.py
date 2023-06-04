import random
import sys
import time
import pandas as pd
import numpy
import argparse
from scoop import futures
from deap import base, creator, tools, algorithms
from LO_properties import sys

from cosimulate import cosimulate, check

def evaluate(individual):
    mu, k, limit = individual
    mu *= 2.0
    individual = mu, k, limit
    return (cosimulate(individual),)


def de_mutation(base_ind_pop, xr1_pop, xr2_pop, F=1):
    l = len(base_ind_pop[0])
    for base_ind, xr1, xr2 in zip(base_ind_pop, xr1_pop, xr2_pop):
        for i in range(l):
            base_ind[i] = base_ind[i] + F * (xr1[i] - xr2[i])
    return base_ind_pop


def binomial_crossover(target_pop, mutant_pop, CR=0.5):
    l = len(target_pop[0])
    for target, mutant in zip(target_pop, mutant_pop):
        k = random.randint(0, l - 1)
        for i in range(l):
            if not (random.random() < CR or i == k):
                mutant[i] = target[i]
    return mutant_pop


def checkBounds(min, max):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in range(len(child)):
                    if child[i] > max:
                        child[i] = max
                    elif child[i] < min:
                        child[i] = min
            return offspring

        return wrapper

    return decorator


def setup(args):
    # Genetic operators
    if args.ga:
        toolbox.register("select", tools.selTournament, tournsize = 2)
        toolbox.register("mate", tools.cxBlend, alpha = 0.5)
        toolbox.decorate("mate", checkBounds(0, 1))
        toolbox.register("mutate", tools.mutGaussian, mu = 0.0, sigma = 0.2, indpb = 0.2)
        toolbox.decorate("mutate", checkBounds(0, 1))
        algorithm = ga
    elif args.de:
        toolbox.register("select", tools.selRandom, k = 3 * 20)
        toolbox.register("mutate", de_mutation, F = 1.0)
        toolbox.decorate("mutate", checkBounds(0, 1))
        toolbox.register("mate", binomial_crossover, CR = 1.0 / 3)
        toolbox.decorate("mate", checkBounds(0, 1))
        algorithm = de
    else:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # To parallelize the evaluation using Scoop
    toolbox.register("map", futures.map)

    # Random seed
    random.seed(int(time.time()))

    # Statistics
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", numpy.mean)
    stats.register("std", numpy.std)
    stats.register("min", numpy.min)
    stats.register("max", numpy.max)

    # Logs
    logbook = tools.Logbook()
    logbook.header = ["gen", "nevals"] + stats.fields

    # Population
    population = toolbox.population(n = 20)

    # Solution
    best = tools.HallOfFame(maxsize = 1)

    return algorithm, population, toolbox, stats, logbook, best


def ga(population, toolbox, stats, logbook, best, n_gen):
    # Setup algorithm
    print('Using genetic algorithm')
    
    # n_gen = 1
    n_best = 1
    cxpb = 1.0
    mutpb = 1.0 / 3

    # Evaluate initial population
    fitnesses = toolbox.map(toolbox.evaluate, population)
    for ind, fit in zip(population, fitnesses):
        ind.fitness.values = fit

    # Update best individual
    best.update(population)

    # Update logs
    record = stats.compile(population)
    logbook.record(gen=0, nevals=len(population), **record)
    print(logbook.stream)

    # Begin the generational process
    for gen in range(1, n_gen):
        # Selection
        offspring = toolbox.select(population, len(population) - n_best)

        # Reproduction
        offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate new individuals
        new_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, new_ind)
        for ind, fit in zip(new_ind, fitnesses):
            ind.fitness.values = fit

        # Keep the best individual of the previous generation
        offspring.extend(best.items)

        # Replace the population
        population[:] = offspring

        # Update the best individual
        best.update(population)

        # Update logs
        record = stats.compile(population)
        logbook.record(gen=gen, nevals=len(new_ind), **record)
        print(logbook.stream)

    # Print best individual
    mu, k, limit = best[0]
    mu *= 2.0
    print("Solución: mu =", mu, " k =", k, " limit =", limit)
    _, _, _, _, _, _, _, _, _, _, cosim_E, _, _, _, indicador, energy, _ = check((mu, k, limit))
    print("Error en la solución:", numpy.linalg.norm(cosim_E - energy) / numpy.sqrt(10001), "J")
    print('Valor del indicador:', indicador[-2], 'J')


def de(population, toolbox, stats, logbook, best, n_gen):
    # Setup algorithm
    # n_gen = 100
    print('Using differential evolutive algorithm')

    # Evaluate initial population
    fitnesses = toolbox.map(toolbox.evaluate, population)
    for ind, fit in zip(population, fitnesses):
        ind.fitness.values = fit

    # Update best individual
    best.update(population)

    # Update logs
    record = stats.compile(population)
    logbook.record(gen=0, nevals=len(population), **record)
    print(logbook.stream)

    # Begin the generational process
    for gen in range(1, n_gen):
        # Selection
        selected = toolbox.select(population)
        l = len(selected)
        base_ind = selected[0 : l // 3]
        xr1 = selected[l // 3 : 2 * l // 3]
        xr2 = selected[2 * l // 3 : l]

        # Mutation
        base_ind = [toolbox.clone(ind) for ind in base_ind]
        mutant = toolbox.mutate(base_ind, xr1, xr2)

        # Crossover
        trial = toolbox.mate(population, mutant)

        # Evaluate new individuals
        fitnesses = toolbox.map(toolbox.evaluate, trial)
        for ind, fit in zip(trial, fitnesses):
            ind.fitness.values = fit

        # Replace the population
        l = len(population)
        for i in range(l):
            if trial[i].fitness > population[i].fitness:
                population[i] = trial[i]

        # Update the best individual
        best.update(population)

        # Update logs
        record = stats.compile(population)
        logbook.record(gen=gen, nevals=len(trial), **record)
        print(logbook.stream)

    # Print best individual
    mu, k, limit = best[0]
    mu *= 2.0
    print("Solución: mu =", mu, " k =", k, " limit =", limit)
    _, _, _, _, _, _, _, _, _, _, cosim_E, _, _, _, indicador, energy, _ = check((mu, k, limit))
    print("Error en la solución:", numpy.linalg.norm(cosim_E - energy) / numpy.sqrt(10001), "J")
    print('Valor del indicador:', indicador[-2], 'J')


# This needs to be here (instead of being in setup()) in order to parallelize with Scoop because it uses Pickle.

# Genes, chromosomes, and population
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)
toolbox = base.Toolbox()
toolbox.register("attr_float", random.random)
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_float, 3)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)
toolbox.register("evaluate", evaluate)

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ga",
        action="store_true",
        help="Genetic algorithm with tournament selection, BLX-\u03B1 crossover and Gaussian mutation",
    )
    parser.add_argument("--de", action="store_true", help="Differential Evolution")
    args = parser.parse_args()

    
    logbooks = list()
    n_gen = 10
    iterations = 10
    p_mu = 0
    p_k = 0
    p_lim = 0
    error = 1e12


    for iteration in range(0, iterations):

        print('Iteration number', iteration + 1, 'out of', iterations)
    
        # Init execution
        algorithm, population, toolbox, stats, logbook, best = setup(args)
        
        # Run algorithm
        algorithm(population, toolbox, stats, logbook, best, n_gen)

        # Copy data
        logbooks.append(logbook)
        p_mu += 1 * best[0][0] / iterations
        p_k += best[0][1] / iterations
        p_lim += best[0][2] / iterations

        # Mejor solución
        _, _, _, _, _, _, _, _, _, _, cosim_E, _, _, _, energy, _ = check((best[0][0], best[0][1], best[0][2]))
        if error > numpy.linalg.norm(cosim_E - energy) / numpy.sqrt(10001):
            
            error = numpy.linalg.norm(cosim_E - energy) / numpy.sqrt(10001)
            p_mu = best[0][0]
            p_k = best[0][1]
            p_lim = best[0][2]

    # Solución media
    print('')
    print('=======================================================================================')
    print("Solución: mu =", p_mu, " k =", p_k, " limit =", p_lim)
    _, _, _, _, _, _, _, _, _, _, cosim_E, _, _, _, indicador, energy, _ = check((p_mu, p_k, p_lim))
    print("Error en la solución:", numpy.linalg.norm(cosim_E - energy) / numpy.sqrt(10001), "J")
    print('Valor del indicador:', indicador[-2], 'J')
    print('=======================================================================================')
    print('')
    
    # Calculate the average
    gn = range(1, 1+n_gen)
    media = numpy.zeros(n_gen)
    minimo = numpy.zeros(n_gen)

    for entry in logbooks:

        _, avg, min = entry.select("gen", "avg", "min")
        
        for generacion in range(0, n_gen):
            media[generacion] += avg[generacion] / iterations
            minimo[generacion] += min[generacion] / iterations

    # Export to csv
    output = pd.DataFrame()
    output['gen'] = gn
    output['media'] = media
    output['min'] = minimo


    output.to_csv('Case1_quality_1_gen.csv', sep = ',', index = False, header = True)




