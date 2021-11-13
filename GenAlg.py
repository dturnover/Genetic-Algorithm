import random


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
    return populationAvg, fittest


# returns two new genomes produced by crossing over the given
# genomes at a random crossover point
def crossover(g1, g2):
    # set the crossover point
    crossoverPoint = random.randrange(0, len(g1))

    # split the first genome at the crossover point
    g1s1 = g1[0:crossoverPoint]
    g1s2 = g1[crossoverPoint:len(g1)]

    # split the second genome at the crossover point
    g2s1 = g2[0:crossoverPoint]
    g2s2 = g2[crossoverPoint:len(g2)]

    # combine the split segments to form offspring
    child1 = g1s1 + g2s2
    child2 = g2s1 + g1s2

    # return the pair of children
    return child1, child2


# returns a new mutated version of the given genome
def mutate(g, mr):
    # access each bit in the genome
    oldg = g
    for i in range(len(g)):
        # determine whether mutation will occur based on mutation rate
        if random.randrange(0, 10000) < mr * 10000:
            # flip the bit
            if g[i] == '0':
                gl = list(g)
                gl[i] = '1'
                g = ''.join(gl)
            elif g[i] == '1':
                gl = list(g)
                gl[i] = '0'
                g = ''.join(gl)

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

    # find the spin's corresponding genome
    for i in range(len(proportions)):
        # if the spin lands on the ith value of proportion then it landed on the (i-1)th value
        if spin >= proportions[i]:
            return i - 1

    # the first index was selected
    return 0


# selects and returns two genomes from the given population using
# fitness-proportionate selection
def selectPair(p, fv):
    # index of the first selection
    location1 = selectSingle(fv, p)
    s1 = p[location1]

    # update structures to remove selected item no. 1 from consideration
    save = [fv[location1], p[location1]]
    del fv[location1], p[location1]

    # index of the second selection
    location2 = selectSingle(fv, p)
    s2 = p[location2]

    # return the deleted values to their respective lists
    fv.insert(location1, save[0])
    p.insert(location1, save[1])

    # return selected pair
    return s1, s2


# returns the new population after performing replacement
def replace(p, fv, g1, g2, x, itr):
    # genome fitness variables
    fg1 = fitness(g1)
    fg2 = fitness(g2)

    # if the selected genomes were crossed
    if x:
        # assign the selected genomes to the two lowest performing genomes
        p[itr] = g1
        p[itr + 1] = g2
        fv[itr] = fg1
        fv[itr + 1] = fg2

        # increase and return the iterator by two, the amount of genomes replaced
        return itr + 2
    else:
        # assign the top selected genome to the lowest performing genome
        if fg1 >= fg2:
            p[itr] = g1
            fv[itr] = fg1
        else:
            p[itr] = g2
            fv[itr] = fg2

        # increase and return the iterator by one, the amount of genomes replaced
        return itr + 1


# is the main GA program, which takes the population size,
# crossover rate, and mutation rate as parameters. The optional
# logFiles parameter is a string specficying the name of a test
# file in which to store the data generated by the GA, for plotting
# purposes. When the GA terminates, this function should return the
# generation at which the string of all ones was found
def runGA(pSize, cRate, mRate, logFile):
    # set the number of generations and create the initial population
    itrFound = 0
    lock = False
    bitLength = 20
    numGenerations = 50
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
        outputStr1 = "Generation: " + str(i + 1) + " | Average fitness: " + str(
            averageFitness) + " | Fittest genome: " + str(fittestGenome)
        outputStr1nl = outputStr1 + '\n'
        print(outputStr1)
        if len(logFile) > 0:
            file.write(outputStr1nl)

        # check if the fittest genome became the optimal string
        if fitness(fittestGenome) == bitLength and lock == False:
            itrFound = i
            lock = True

        # half the population size iterations since two genomes selected at a time
        itr = 0
        for j in range(int(pSize) // 2):
            # select two genomes from the population based on fitness-proportionate selection
            genome1, genome2 = selectPair(population, fitnessValues)

            # determine whether crossover will occur based on crossover rate
            crossed = False
            if random.randrange(0, 100) < cRate * 100:
                genome1, genome2 = crossover(genome1, genome2)
                crossed = True

            # mutate the selected genomes
            genome1 = mutate(genome1, mRate)
            genome2 = mutate(genome2, mRate)

            # replace the worst genome(s) with the selected genome(s)
            itr = replace(population, fitnessValues, genome1, genome2, crossed, itr)

    if itrFound > 0:
        outputStr2 = "Optimal solution found in generation " + str(itrFound + 1) + '\n'
        print(outputStr2)
        if len(logFile) > 0:
            file.write(outputStr2)

    if len(logFile) > 0:
        file.close()


fileName = input("Enter a log file name or press 'enter' to skip: ")
runGA(100, 0.7, 0.001, fileName)
