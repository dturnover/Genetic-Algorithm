import random
import matplotlib.pyplot as plt
import numpy
import numpy as np

# Debug mode, graph mode and print mode toggles
Debug = False
Graph = False
Print = True


# selection sort algorithm
def selectionSort(lscores, lreg):
    for i in range(len(lscores)):
        # Find the minimum element in unsorted array
        min_idx = i
        for j in range(i + 1, len(lscores)):
            if lscores[min_idx] > lscores[j]:
                min_idx = j

        # Swap the found minimum element with the first element
        lscores[i], lscores[min_idx] = lscores[min_idx], lscores[i]
        lreg[i], lreg[min_idx] = lreg[min_idx], lreg[i]

    # return the reordered fitness scores and population arrays
    return lscores, lreg


# returns a random individual(genome) (bit string) of a given length
def randomGenome(l):
    bitString = ''

    # for-loop that generates a bitsring using random bits
    for i in range(l):
        # generate random bit
        randomBit = str(random.randrange(0, 2))
        bitString += randomBit

    # return bit string
    return bitString


# returns a new randomly created population of the specified size,
# represented as a list ofgenomes of the specified length
def makePopulation(s, l):
    # empty population definition
    population = []

    # for-loop that fills the population array with a bitstring
    for i in range(s):
        bitString = randomGenome(l)
        population.append(bitString)

    # return the newly generated population
    return population


# returns a list of corresponding fitness values for a population
def populationFitness(p):
    fitnessValues = []
    for genome in p:
        # calculate the fitness value of the genome
        f = fitness(genome)

        # append it to the new list
        fitnessValues.append(f)

    # return the population's fitness values
    return fitnessValues


# returns the fitness value of a genome
def fitness(g):
    # loop through the bitsring
    count = 0
    for c in g:
        # count the number of ones
        if c == '1':
            count += 1

    # return the count
    return count


# returns a pair of values: the average fitness of the population as
# a whole and the fitness of the best individual in the population
def evaluateFitness(p):
    # variables to keep track of the fittest genome and population avg
    fittest = ''
    fHigh = 0
    fTotal = 0

    # access each member of the population
    for genome in p:
        # calculate fitness of genome
        f = fitness(genome)
        # check if it is the fittest thus far
        if f > fHigh:
            fittest = genome
            fHigh = f
        # update population's fitness score
        fTotal += f
    populationAvg = fTotal / len(p)

    # return average fitness and fittest genome
    return round(populationAvg, 1), fittest


# returns two new genomes produced by crossing over the given
# genomes at a random crossover point
def crossover(g1, g2):
    # set the crossover point
    crossoverPoint = random.randrange(0, len(g1))

    # debug statement
    if Debug:
        print("Performing crossover:\n- Crossover point:", crossoverPoint)
        print("- Parent genomes:", g1, ",", g2)

    # split the first genome at the crossover point
    g1s1 = g1[0:crossoverPoint]
    g1s2 = g1[crossoverPoint:len(g1)]

    # split the second genome at the crossover point
    g2s1 = g2[0:crossoverPoint]
    g2s2 = g2[crossoverPoint:len(g2)]

    # combine the split segments to form offspring
    child1 = g1s1 + g2s2
    child2 = g2s1 + g1s2

    # debug statement
    if Debug:
        print("- Child one:", child1)
        print("- Child two:", child2)

    # return the pair of children
    return child1, child2


# returns a new mutated version of the given genome
def mutate(g, mr):
    # access each bit in the genome
    oldg = g
    for i in range(len(g)):
        # determine whether mutation will occur based on mutation rate
        if random.randrange(0, 10000) < mr * 10000:
            # debug statement
            if Debug:
                print("------------------------------------------------\nMutating:\n- Before:", g)
            # flip the bit
            if g[i] == '0':
                gl = list(g)
                gl[i] = '1'
                g = ''.join(gl)
            elif g[i] == '1':
                gl = list(g)
                gl[i] = '0'
                g = ''.join(gl)

            # debug statement
            if Debug:
                print("- After: ", g)

    # return the genome
    return g


# selects and returns one of the two gnomes from the population
# using fitness-proportionate selection
def selectSingle(fv, p):
    # fill the proportions array with dividen values
    proportions = []
    fTotal = 0
    for value in fv:
        fTotal += value
        proportions.append(fTotal)

    # obtain random number ('spin' the wheel)
    spin = random.randrange(0, fTotal)

    # debug statement
    if Debug:
        print("- Test spin:", spin)

    # find the spin's corresponding genome
    for i in range(len(proportions)):
        # if the spin lands on the ith value of proportion then it landed on the (i-1)th value
        if spin <= proportions[i]:
            # debug statement
            if Debug:
                print("- Test proportions:", proportions[i], "i:", i)

            return i - 1

    # the first index was selected
    return 0


# selects and returns two genomes from the given population using
# fitness-proportionate selection
def selectPair(p, fv):
    # debug statement
    if Debug:
        print("\nSelecting parents:")

    # index of the first selection
    location1 = selectSingle(fv, p)
    s1 = p[location1]

    # debug statement
    if Debug:
        print("- Test locations 1:", location1, "s1:", s1)

    # update structures to remove selected item no. 1 from consideration
    save = [fv[location1], p[location1]]
    del fv[location1], p[location1]

    # index of the second selection
    location2 = selectSingle(fv, p)
    s2 = p[location2]

    # debug statement
    if Debug:
        print("- Test locations 2:", location2, "s2:", s2)

    # return the deleted values to their respective lists
    fv.insert(location1, save[0])
    p.insert(location1, save[1])

    # return selected pair
    return s1, s2


# returns the new population after performing replacement
def replace(p, fv, g1, g2, x):
    # genome fitness variables
    fg1 = fitness(g1)
    fg2 = fitness(g2)

    # if the selected genomes were crossed
    if x:
        # debug statement
        if Debug:
            print("Replacing:", p[0], "and", p[1], " with", g1, "and", g2, ", respective fitnesses:", fitness(p[0]), fitness(p[1]), "->", fg1, fg2)

        # assign the selected genomes to the two lowest performing genomes
        p[0] = g1
        p[1] = g2
        fv[0] = fg1
        fv[1] = fg2
    else:
        # assign the top selected genome to the lowest performing genome
        if fg1 >= fg2:
            # debug statement
            if Debug:
                print("Replacing:", p[0], "with", g1, ", respective fitnesses:", fitness(p[0]), "t->",  fg1)

            p[0] = g1
            fv[0] = fg1
        else:
            # debug statement
            if Debug:
                print("Replacing:", p[0], "with", g2, ", respective fitnesses:", fitness(p[0]), "->",  fg2)

            p[0] = g2
            fv[0] = fg2


# is the main GA program, which takes the population size,
# crossover rate, and mutation rate as parameters. The optional
# logFiles parameter is a string specficying the name of a test
# file in which to store the data generated by the GA, for plotting
# purposes. When the GA terminates, this function should return the
# generation at which the string of all ones was found
def runGA(pSize, cRate, mRate, logFile):
    # set the number of generations and create the initial population
    itrFound = 0
    bitLength = 20
    numGenerations = 1000
    population = makePopulation(pSize, bitLength)

    # open a file to write to
    file = 0
    if len(logFile) > 0:
        file = open(logFile, 'w')

    # each iteration simulates a single generation
    for i in range(numGenerations):
        # sort the population and fitnessValues arrays in ascending order
        fitnessValues = populationFitness(population)
        selectionSort(fitnessValues, population)

        # report the fitness of the current generation
        averageFitness, fittestGenome = evaluateFitness(population)
        outputStr1 = "Generation: " + str(i + 1) + " | Average fitness: " + str(averageFitness) + " | Fittest genome: " + str(fittestGenome)
        outputStr1nl = outputStr1 + '\n'
        if Print:
            print(outputStr1)
        if len(logFile) > 0:
            file.write(outputStr1nl)

        # check if the fittest genome became the optimal string
        if fitness(fittestGenome) == bitLength:
            itrFound = i
            break

        # select two genomes from the population based on fitness-proportionate selection
        genome1, genome2 = selectPair(population, fitnessValues)

        # determine whether crossover will occur based on crossover rate
        crossed = False
        if random.randrange(0, 100) < cRate * 100:
            genome1, genome2 = crossover(genome1, genome2)
            crossed = True
        else:
            # debug statement
            if Debug:
                print("NOT CROSSED")

        # mutate the selected genomes
        genome1 = mutate(genome1, mRate)
        genome2 = mutate(genome2, mRate)

        # replace the worst genome(s) with the selected genome(s)
        replace(population, fitnessValues, genome1, genome2, crossed)

    # this statement checks if the optimal solution was found at any point during execution
    if itrFound > 0 and Print:
        outputStr2 = "Optimal solution found in generation " + str(itrFound + 1) + '\n'
        print(outputStr2)
        if len(logFile) > 0:
            file.write(outputStr2)

    # close the log file if it was opened
    if len(logFile) > 0:
        file.close()

    # return the iteration at which the optimal solution was found
    return itrFound


# graph mode runs GA many times and plots the adjusted variable vs output
if Graph:
    # population size vs number of generations
    numGenerationsPs = [[], []]
    for populationSize in range(10, 60, 5):
        # thirty to make use of central limit theorem
        for i in range(30):
            numGenerationsPs[0].append(populationSize)
            numGenerationsPs[1].append(runGA(populationSize, 0.9, 0.01, ""))
    numGenerationsPsX = np.array(numGenerationsPs[0])
    numGenerationsPsY = np.array(numGenerationsPs[1])

    # crossover rate vs number of generations
    numGenerationsCr = [[], []]
    for crossoverRate in range(10):
        crossoverRateF = float(crossoverRate)/10
        # thirty to make use of central limit theorem
        for j in range(30):
            numGenerationsCr[0].append(crossoverRateF)
            numGenerationsCr[1].append(runGA(30, crossoverRateF, 0.01, ""))
    numGenerationsCrX = np.array(numGenerationsCr[0])
    numGenerationsCrY = np.array(numGenerationsCr[1])

    # mutation rate vs number of generations
    numGenerationsMr = [[], []]
    for mutationRate in range(10):
        mutationRateF = float(mutationRate)/100
        # thirty to make use of central limit theorem
        for k in range(30):
            numGenerationsMr[0].append(mutationRateF)
            numGenerationsMr[1].append(runGA(30, 0.9, mutationRateF, ""))
    numGenerationsMrX = np.array(numGenerationsMr[0])
    numGenerationsMrY = np.array(numGenerationsMr[1])

    # Display scatterplots
    plt.figure(figsize=(7, 7))
    plt.xlabel('Population Size')
    plt.ylabel('Number of Generations')
    plt.scatter(numGenerationsPsX, numGenerationsPsY)
    plt.title('Population Size vs Number of Generations')

    plt.figure(figsize=(7, 7))
    plt.xlabel('Crossover Rate')
    plt.ylabel('Number of Generations')
    plt.scatter(numGenerationsCrX, numGenerationsCrY)
    plt.title('Crossover Rate vs Number of Generations')

    plt.figure(figsize=(7, 7))
    plt.xlabel('Mutation Rate')
    plt.ylabel('Number of Generations')
    plt.scatter(numGenerationsMrX, numGenerationsMrY)
    plt.title('Mutation Rate vs Number of Generations')

    plt.show()
else:
    # run normally, just one time when not in graph mode
    fileName = input("Enter a log file name or press 'enter' to skip: ")
    runGA(30, 0.9, 0.01, fileName)
