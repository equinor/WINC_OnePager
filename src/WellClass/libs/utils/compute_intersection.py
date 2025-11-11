import numpy as np

def compute_intersection(x: np.ndarray, y1:np.ndarray, y2: np.ndarray):
    
    x = np.asarray(x, dtype=float)
    y1 = np.asarray(y1, dtype=float)
    y2 = np.asarray(y2, dtype=float)

    # Clean data
    m = np.isfinite(x) & np.isfinite(y1) & np.isfinite(y2)
    x, y1, y2 = x[m], y1[m], y2[m]


    # Ensure x is sorted (required for interpolation to make sense)
    s = np.argsort(x)
    x, y1, y2 = x[s], y1[s], y2[s]


    # Difference between curves
    d = y1 - y2

    # Indices where sign changes (crossing zero)
    idx = np.where(np.diff(np.sign(d)))[0]

    # Linear interpolation for more accurate intersection points
    xi = x[idx] - d[idx] * (x[idx+1] - x[idx]) / (d[idx+1] - d[idx])
    yi = y1[idx] - d[idx] * (y1[idx+1] - y1[idx]) / (d[idx+1] - d[idx])

    if xi.shape[0] == 0:
        print('Lines do not intersect')
        intersect_x = np.nan
        intersect_y = np.nan

    elif xi.shape[0] > 1:
        intersect_x = xi[np.argmax(xi)]
        intersect_y = yi[np.argmax(xi)]

    else:
        intersect_x = xi[0]
        intersect_y = yi[0]
    
    
    return intersect_x, intersect_y
