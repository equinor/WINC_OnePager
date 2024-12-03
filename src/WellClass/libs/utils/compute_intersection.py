import numpy as np

def compute_intersection(x: np.ndarray, y1:np.ndarray, y2: np.ndarray):

        """
        Function that computes intersection between two curves.
        takes 4 parameters:
        x = depth vectors. By defaultf they share the same sampling. 
        y1, y2 = Sh_min and fluid pressure vectors.

        """
        #Convert to float
        x = x.astype(float)
        y1 = y1.astype(float)
        y2 = y2.astype(float)
        
        # Filter nan values in vector
        
        y_filter = ~np.isnan(y2)
        x =   x[y_filter]
        y1 = y1[y_filter]
        y2 = y2[y_filter]

        #Retrieve closest index
        idx = np.argwhere(np.diff(np.sign(y1 - y2))).flatten()

        #Define subsets with x, y1 and y2 values at retrieved idx
        x_subset = [x[idx]]
        y1_subset = [y1[idx]]
        y2_subset = [y2[idx]]

        #offsets to evaluate prior or forward points in vector
        offsets = [-1, 1]

        #iterate over offsets and check if intervals intersect.
        for offset in offsets:
                y10 = y1[idx]
                y11 = y1[idx + offset]
                y20 = y2[idx]
                y21 = y2[idx + offset]

                #check if pressure (y1) at idx is higher than pressure (y2) at idx
                verif_0 = y10>y20

                #Verify intersection
                if verif_0:
                        verif_1 = y11<y21 

                else:
                        verif_0 = y10<y20
                        verif_1 = y11>y21

                #Find crossing indexes
                if verif_0 and verif_1:
                        x_subset.append(x[idx+offset])
                        y1_subset.append(y11)
                        y2_subset.append(y21)


        #Compute linear expression for each subset
        #compute slopes m1 and m2
        if len(y2_subset)>1 and len(y1_subset)>1:
                m1 = (y1_subset[0]-y1_subset[1])/(x_subset[0]-x_subset[1])
                m2 = (y2_subset[0]-y2_subset[1])/(x_subset[0]-x_subset[1])

                #compute intercepts b1 and b2
                b1 = y1_subset[0]-x_subset[0]*m1
                b2 = y2_subset[0]-x_subset[0]*m2

                #find intersection
                intersect_x = (-(b1-b2)/(m1-m2))[0]
                intersect_y = (m1*intersect_x+b1)[0]

        else:
                print('Lines do not intersect')
                intersect_x = np.nan
                intersect_y = np.nan

        return intersect_x, intersect_y
