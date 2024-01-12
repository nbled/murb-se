import subprocess

optims = [
#    "cpu+naive",
    "cpu+optim",
    "cpu+omp",
    "cpu+simd",
    "cpu+simd+omp",
    "cpu+simd+pthread",
    "cpu+barnesHut",
    "cpu+barnesHut+omp",
    "gpu"
]

def run_murb(optim, nbodies, iter):
    max_fps = 0
    min_ms = 0

    for _ in range(10):

        proc = subprocess.Popen([
            "./bin/murb", "--im", optim, 
            "-n", str(nbodies), 
            "-i", str(iter), 
            "--nv"], 
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, env={"OMP_NUM_THREADS": "6"})
        output = proc.stdout.read().decode().splitlines()
        status = output[-1].split()
        
        ms, fps = float(status[-4]), float(status[-2][1:])
        if fps > max_fps:
            max_fps = fps
            min_ms = ms

        if ms > 10000:
            break

    return min_ms, max_fps

nbodies_list = [500, 1000, 2000, 3000, 4000, 5000]

"""
csv = open("fps_per_bodies.csv", "wt")
csv.write("optim,nbodies,ms,fps\n")

for nbodies in nbodies_list:
    for optim in optims:
        ms, fps = run_murb(optim, nbodies, 1000)
        print(f"{optim},{nbodies},{ms},{fps}")
        csv.write(f"{optim},{nbodies},{ms},{fps}\n")

csv.close()

 """


csv = open("../fps_per_theta.csv", "at")
csv.write("optim,nbodies,ms,fps\n")



optim = "cpu+barnesHut"
theta = 1
for nbodies in nbodies_list:
    ms, fps = run_murb(optim, nbodies, 30)
    print(f"B-H {theta},{nbodies},{ms},{fps}")
    csv.write(f"B-H {theta},{nbodies},{ms},{fps}\n")


