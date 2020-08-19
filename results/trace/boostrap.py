import numpy as np 

def CV(mean, x):
    return np.var(x, ddof=1)/mean**2


A1 = np.genfromtxt("ss_mRNA__10_10_0.01_0.01_1.41_1.dat")
A05 = np.genfromtxt("ss_mRNA__10_10_0.01_0.01_0.14_0.05.dat")

a1_mean = []
a1_cv= []
a05_mean= []
a05_cv= []

for i in range(100):
    a1 = np.random.choice(A1, len(A1), replace=True)
    a1_mean.append(np.mean(a1)) 
    a1_cv.append(CV(a1_mean[-1], a1))

    a05 = np.random.choice(A05, len(A1), replace=True)
    a05_mean.append(np.mean(a05)) 
    a05_cv.append(CV(a05_mean[-1], a05))

print("A1 mean", np.mean(a1_mean), np.std(a1_mean), )
print("A1 CV2", np.mean(a1_cv), np.std(a1_cv), )

print("A05 mean", np.mean(a05_mean), np.std(a05_mean), )
print("A02 CV2", np.mean(a05_cv), np.std(a05_cv), )



# A1 mean 70.17283900000001 0.6317590455062759
# A1 CV2 1.0239243347230012 0.02031471765164097
# A05 mean 71.14697600000001 0.7060320154383931
# A02 CV2 0.7957485888690067 0.012416012993820197


