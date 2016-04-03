import numpy as np
import sys

def main():
  ccmfn = sys.argv[1]
  admfn = sys.argv[2]

  ccm = np.loadtxt(ccmfn)
  adm = np.loadtxt(admfn)

  assert np.array_equal(ccm, ccm.T)

  summed = ccm + adm + adm.T
  assert np.sum(summed[summed < 0.0]) == 0
  more_than_one = summed[summed > 1.0]
  # Floating-point addition errors may result in some elements in summed being
  # slightly more than one, so can't do same assert we did for < 0.
  assert np.allclose(more_than_one, 1.0)

main()
