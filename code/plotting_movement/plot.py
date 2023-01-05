import numpy as np
import matplotlib.pyplot as plt

DATA = np.loadtxt("../n_bodies/output_pthread2.txt", delimiter=",")
ranges = [np.min(DATA[:, 2:5], axis=0), np.max(DATA[:, 2:5], axis=0)]
N = (DATA[:,0] == 0).sum()
plt.figure(figsize=(8,4))
# ax = plt.axes(projection='3d')

for i in range(len(DATA)//N):
    # col = DATA[N*i:N*(i+1), 4]

    if i % 100 == 0:
        col = np.zeros((N, 4))
        col[:, 3] = 1
        col[:, 2] = ((DATA[N*i:N*(i+1), 4] - ranges[0][2]) / ranges[1][2]) 

        plt.cla()
        plt.title(f"Iteration: {i}")

        plt.xlim([ranges[0][0],ranges[1][0]])
        plt.ylim([ranges[0][1], ranges[1][1]])
        plt.scatter(DATA[:N*(i+1), 2],DATA[:N*(i+1), 3], s=1, color=[.7, .7, 1])

        plt.scatter(DATA[N*i:N*(i+1), 2], DATA[N*i:N*(i+1), 3], s=20, c=col)
        plt.pause(0.00001)
plt.show()
