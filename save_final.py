from decimal import *
from math import comb

getcontext().prec = 400

m = 8000000000
k = 120
num_iterations = 400


def eval_derivative_approx(m, k, p):
    first = Decimal(p) ** (k)
    if k == m:
        first -= ((1 - Decimal(p)) ** (m - 1)) / (Decimal(p) ** (m - k - 1))
    for n in range(1, min(m, k + 1)):
        first -= (
            Decimal(comb(m, n))
            * (Decimal(1 - p) ** (n - 1))
            * (Decimal(p) ** (k - n + 1))
        )
        first += (
            Decimal(comb(m, n))
            / n
            * (m - n)
            * (Decimal(1 - p) ** (n))
            * (Decimal(p) ** (k - n))
        )
    return first


def eval_derivative(m, p):
    first = Decimal(p) ** (m - 1)
    first -= (1 - Decimal(p)) ** (m - 1)
    for n in range(1, m):
        first -= (
            Decimal(comb(m, n)) * (Decimal(1 - p) ** (n - 1)) * (Decimal(p) ** (m - n))
        )
        first += (
            Decimal(comb(m, n))
            / n
            * (m - n)
            * (Decimal(1 - p) ** (n))
            * (Decimal(p) ** (m - n - 1))
        )
    return first


# Assume m >> k
def eval_actual_approx(m, k, p):
    first = (Decimal(p) ** (m)) / m
    for n in range(1, k + 1):
        first += (
            Decimal(comb(m, n)) / n * (Decimal(1 - p) ** (n)) * (Decimal(p) ** (m - n))
        )
    return first


upper = 1 - (Decimal(1) / m)
assert eval_derivative_approx(m, k, upper) < 0
lower = upper
while eval_derivative_approx(m, k, lower) < 0:
    lower -= Decimal(1) / m
    assert lower > 0

for _ in range(num_iterations):
    mid = (lower + upper) / 2
    result = eval_derivative_approx(m, k, mid)
    if result < 0:
        upper = mid
    else:
        lower = mid

print(lower, eval_actual_approx(m, k, lower))
