#!/usr/bin/python3
# -*- coding: utf-8 -*-

from math import log10
from maths.primes import primes_range
import matplotlib.pyplot as plt


def inverse_primorial_graph(bound=100):
    """Draws the graphic of the inverse primorial,
    ie the prime bound p = f(n) such that #p = log2(n)"""
    primes = primes_range(2, bound)

    l = len(primes)
    x = [0] * (l + 1)
    for i in range(l):
        x[i + 1] = log10(primes[i]) + x[i]

    majoration = [9 + 2.6 * k for k in x]

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.set_title("Greatest prime p needed such that #(p-1) < n ≤ #p")
    ax.set_xlabel("log(n)")
    ax.set_ylabel("greatest prime p needed")

    ax.set_xlim(left=0, right=int(x[-1]))
    ax.set_ylim(bottom=0, top=bound)

    plt.grid()

    primes_step = []
    for e in primes:
        primes_step.append(e)
        primes_step.append(e)

    x_steps = []
    for e in x:
        x_steps.append(e)
        x_steps.append(e)
    x_steps.pop()

    ax.plot(x_steps, [0] + primes_step, label="Inverse primorial")
    ax.plot(x, majoration, label=r"Majoration : 2.6 * log(n) + 9")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    inverse_primorial_graph()