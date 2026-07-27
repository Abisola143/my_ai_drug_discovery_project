from scipy.special import comb
def calculateProbablility(k, m, n):
    # finding the total population
    totalPop = k + m + n
    # finding all the potential combos of pairings
    totalCombos = comb(totalPop, 2)
    # finds all the valid combos that could make dominant alleles
    validCombos = comb(k,2) + k*m + k*n + .5*m*n + .75*comb(m, 2)

    probability = validCombos/totalCombos
    print(probability)


with open(r"C:\Users\Akinyemi\Downloads\rosalind_iprb (1).txt", "r") as file:
    line = file.readline().split()
    # iterating through the list of words generated
    k,m,n = [int(n) for n in line]
    print(k,m,n)
file.close()

calculateProbablility(k, m, n)