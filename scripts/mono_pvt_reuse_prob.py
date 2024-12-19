import math
import matplotlib.pyplot as plt

K = 126 # 126, 200, 274
N = 2*K

def n_choose_k(n,k):
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

# Full space of first half of the permutation is N choose K
full_space = n_choose_k(N,K)

probs = []
for p in range(K):
    R = p

    ####################################
    ### Odds of exactly R reusable pivots
    ####################################
    ## Choice from left side is K choose R
    space = n_choose_k(K,R)

    ## Choice from right side is K choose (K-R)
    space *= n_choose_k(K,K-R)

    probs.append(space/full_space)

print("Probability of exactly R reusable pivots")
print(probs)
print("----------------")

prob_cum = []
s = 0
for i in range(K):
    s += probs[i]
    prob_cum.append(s)

print("Cumulative probability of exactly R reusable pivots")
print(prob_cum)
print("----------------")

print("Cumulative probability of exactly R reusable pivots (as powers of 2)")
print([math.log2(x) for x in prob_cum])
print("----------------")

for i in range(K//2):
    print(f"R={i}, Prob of fewer than R pivots=2^{math.log2(prob_cum[i])}")


fig = plt.figure()
plt.plot(probs)
plt.title("Probablility R Reusable Pivots")
plt.xlabel("R")
fig.savefig(f'prob_{N}_{K}.png', dpi=fig.dpi)

fig = plt.figure()
plt.plot(prob_cum)
plt.title("Cumulative Probablility R Reusable Pivots")
plt.xlabel("R")
fig.savefig(f'prob_cum_{N}_{K}.png', dpi=fig.dpi)