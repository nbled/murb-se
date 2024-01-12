import matplotlib.pyplot as plt

with open("fps_per_theta.csv", "rt") as fp:
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


# with open("fps_per_optim_30000.csv", "rt") as fp:
    # data = fp.read()
# 
# lines = data.splitlines()[1:]
# 
# curves = {}
# for line in lines:
    # optim, nbodies, ms, fps = line.split(",")
    # if optim.startswith("cpu"):
        # optim = optim[4:]
#         
    # curves[optim] = round(float(fps),2)
# 
# 
# plt.xlabel("optim")
# #plt.yscale("log")
# plt.ylabel("FPS")
# 
# 
# bars = plt.bar(list(range(len(curves))),list(curves.values()))
# for i, rect in enumerate(bars):
    # height = rect.get_height()
    # plt.text(rect.get_x() + rect.get_width() / 2.0, height, list(curves.keys())[i], ha='center', va='bottom')
# 
# plt.legend()
# plt.show()
