import matplotlib.pyplot as plt

with open("fps_per_bodies.csv", "rt") as fp:
    data = fp.read()

lines = data.splitlines()[1:]

curves = {}
for line in lines:
    optim, nbodies, ms, fps = line.split(",")
    if optim not in curves:
        curves[optim] = []
    curves[optim].append(float(fps))

nbodies_list = [500, 1000, 2000, 3000, 4000, 5000]

plt.xlabel("Number of bodies")
plt.yscale("log")
plt.ylabel("FPS")

for optim, points in curves.items():
    plt.plot(nbodies_list, points, label=optim)

plt.legend()
plt.show()
