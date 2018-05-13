#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 16:12:26 2017

@author: caoxiya
"""

import numpy as np
import matplotlib.pyplot as plt
import random

def readResult(filename):
    file = open(filename)
    results = []
    result = []
    line = file.readline()
    while 1:
        while 1:
            line = file.readline()
            try:
                result.append(float(line))
            except:
                results.append(result)
                result = []
                print line
                break
            if not line:
                break
        if not line:
            break
        pass # do something
    return results

Aresults = readResult("./dataset/result/p08_annel_result.txt")
results = readResult("./dataset/result/p08_ga_result.txt")
Rresults = readResult("./dataset/result/p08_rand_result.txt")
HCresults = readResult("./dataset/result/p08_hc_result.txt")


HCresult = []
result = []
Rresult = []
Aresult = []
for j in range(0,3):
    hc = []
    ga = []
    rand = []
    anneal = []
    for i in range(0, 100):
        hc.append(HCresults[j][i])
        ga.append(results[j][i * 25])
        rand.append(Rresults[j][i * 25])
        anneal.append(Aresults[j][i * 25])
    HCresult.append(hc)
    result.append(ga)
    Rresult.append(rand)
    Aresult.append(anneal)



HCresults = np.array(HCresult)
HCresults_mean = sum(HCresults,2)/3
HCresults_std = np.std(HCresults.transpose(),1)
Rresults = np.array(Rresult)
Rresults_mean = sum(Rresults,2)/3
Rresults_std = np.std(Rresults.transpose(),1)
results = np.array(result)
results_mean = sum(results,2)/3
results_std = np.std(results.transpose(),1)
Aresults = np.array(Aresult)
Aresults_mean = sum(Aresults,2)/3
Aresults_std = np.std(Aresults.transpose(),1)

plt.plot(results_mean, color = 'g', label = 'Genetic Algorithm')
plt.plot(Rresults_mean, color = 'b', label = 'Random Search')
plt.plot(HCresults_mean, color = 'r', label = 'Hill Climber')
plt.legend(loc='best')
plt.show()
x = np.array(range(0,100))
plt.errorbar(x,results_mean,results_std,color = 'g',label = 'GeneticAlgorithm')
plt.errorbar(x,HCresults_mean, HCresults_std,color = 'r',label = 'HillClimber')
plt.errorbar(x,Rresults_mean, Rresults_std,color = 'b', label = 'RandomSearch')
#plt.errorbar(x,Aresults_mean, Aresults_std,color = 'y', label = 'Simulated Annealing')
plt.grid(True)
plt.legend(loc='best')
plt.show()
