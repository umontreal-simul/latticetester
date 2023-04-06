from random import randrange
primes = [1021, 1048573, 1073741827, 1099511627791, 1125899906842597, 18446744073709551629]
for prime in primes:
    for j in range(5, 101, 5):
        for i in range(10):
            a = []
            for k in range(j):
                a.append([])
                for l in range(j):
                    a[k].append(randrange(prime))
            with open(str(prime) + '_' + str(j) + '_' + str(i) + '.dat', 'w') as f:
                f.write(str(j) + '\n')
                for k in range(j):
                    for l in range(j):
                        f.write(str(a[k][l]) + ' ')
                    f.write('\n')
