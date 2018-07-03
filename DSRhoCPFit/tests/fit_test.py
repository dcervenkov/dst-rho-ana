import subprocess
import sys

print("Running a full fit DSRhoCPFit test")
return_code = subprocess.call(["./DSRhoCPFit", "--efficiency-model=6", "--mixing", "--events=1000",
                               "--fix=apa,a0,ata,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb", "--",
                               "tests/signalMC_data.root", "tests/current_result"],
                               cwd="../.")

reference_result = ""
current_result = ""

with open("reference_result", "r") as f:
    reference_result = f.readline()

with open("current_result", "r") as f:
    current_result = f.readline()

if return_code == 0 and reference_result == current_result:
    print("The result matches the reference result.")
    sys.exit(0)
else:
    print("The result does NOT match the reference result!")
    sys.exit(1)
