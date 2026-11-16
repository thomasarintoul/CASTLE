import numpy as np

bins = np.fromfile('sf_L6_sur1_tng_full_test_bin_edges.bin', dtype=np.float64)
counts = np.fromfile('sf_L6_sur1_tng_full_test_counts.bin', dtype=np.float64)
vsf = np.fromfile('sf_L6_sur1_tng_full_test_velocity_VectorDifference_order_2.bin', dtype=np.float64)
tsf = np.fromfile('sf_L6_sur1_tng_full_test_T_ScalarDifference_order_2.bin', dtype=np.float64)

first_10_vsf = vsf[:10] / counts[:10]
first_10_tsf = tsf[:10] / counts[:10]

expected_vsf = np.array([4.08124262, 3.95935851, 3.89760771, 3.97661365, 4.00049998, 4.00605026, 3.97096193, 3.95534304, 4.04506672, 4.01354214])
expected_tsf = np.array([9.75200909e+08, 1.14324642e+09, 3.02294464e+09, 2.84802940e+09, 2.89427911e+09, 6.67378846e+09, 3.29616119e+09, 3.53919471e+09,  6.64230469e+09, 5.40784152e+09])

print("First 10 values of the velocity structure function (vsf):")
print(first_10_vsf)
diff = np.abs(first_10_vsf - expected_vsf) / expected_vsf
if np.all(diff < 1e-6):
    print("Velocity structure function values match expected values.")
else:
    print("Velocity structure function values do not match expected values.")


print('-----')
print("First 10 values of the temperature structure function (tsf):")
print(first_10_tsf)
diff_tsf = np.abs(first_10_tsf - expected_tsf) / expected_tsf
if np.all(diff_tsf < 1e-6):
    print("Temperature structure function values match expected values.")
else:
    print("Temperature structure function values do not match expected values.")

