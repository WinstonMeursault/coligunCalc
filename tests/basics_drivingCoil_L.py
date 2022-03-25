import sys
sys.path.append(".")

from basics import drivingCoil


test = drivingCoil(rdi=0.043, rde=0.06, ld=0.08, n=77, x0=0)

# 预期输出:≈0.000428
print(test.L)

if abs(test.L - 0.000428) < 0.00002:
    print("[TEST] basics_drivingCoil_L:PASS")
else:
    print("[TEST] basics_drivingCoil_L:FAIL")
