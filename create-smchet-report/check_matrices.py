import numpy as np
import sys

def main():
  ccmfn = sys.argv[1]
  admfn = sys.argv[2]

  ccm = np.loadtxt(ccmfn)
  adm = np.loadtxt(admfn)

  assert np.array_equal(ccm, ccm.T)

  summed = ccm + adm + adm.T
  assert np.sum(summed[summed < 0]) == 0
  assert np.sum(summed[summed > 1]) == 0

main()
