import matplotlib.pyplot as plt

u = []
v = []

"""
0 er plot for 1.c
1 er plot for 1.e
"""
b = 0

if b == 0:
    f = open("verdier.txt", "r")
    n = f.readline()
    t = n.split()[1]
    n = n.split()[0]
    for i in range(int(n)-1):
        n = f.readline()
        n = n.split()
        u.append(float(n[0]))
        v.append(float(n[1]))


elif b == 1:
    f = open("verdier2.txt", "r")
    n = f.readline()
    t = n.split()[1]
    n = n.split()[0]
    for i in range(int(n)):
        a = float(f.readline())
        u.append(a)

    for i in range(int(n)):
        a = float(f.readline())
        v.append(a)
        print(v[i])

else:
    print("Velg enten 1 eller 0 for b")

plt.plot(u, v)
plt.title("t = " + t)
plt.legend()
plt.show()
