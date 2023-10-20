import pytest
import numpy as np

from src.GaP.libs.models import (
    PipeCementModel,
    DepthModel,
    ElemModel,
)

@pytest.fixture(scope='function')
def casing_fixture():
    """ fixture for casing
    """
    # conductor casing info
    cond_casing_geom = PipeCementModel(ID=0.762, 
                                       pipe=DepthModel(strt_depth=312.0, end_depth=371.0, perm=10000), 
                                       oph=DepthModel(strt_depth=4.0, end_depth=376.0, perm=None), 
                                       cement=DepthModel(strt_depth=312.0, end_depth=371.0, perm=5), 
                                       type='conductor')
 
    return cond_casing_geom

@pytest.fixture(scope='function')
def barrier_fixture():
    """ fixture for barrier
    """

    # barrier info
    barrier_geom = ElemModel(ID=0.3397, 
                             pipe=DepthModel(strt_depth=352.0, end_depth=560.0, perm=0.5), 
                             type='barrier')
    
    return barrier_geom

@pytest.fixture(scope='function')
def LGR_fixture():
    """ fixture for LGR grid
    """
    # xy grid sizes
    LGR_sizes_xy = [94.0, 
  5.0,
  0.5,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.05,
  0.5,
  5.0,
  94.0]
    
    # depth grid
    LGR_depths = np.array([  34.8,   65.6,   96.4,  127.2,  158. ,  188.8,  219.6,  250.4,
         281.2,  312. ,  322. ,  332. ,  342. ,  352. ,  362. ,  372. ,
         382. ,  392. ,  402. ,  412. ,  422. ,  432. ,  442. ,  452. ,
         462. ,  472. ,  482. ,  492. ,  502. ,  512. ,  522. ,  532. ,
         542. ,  552. ,  562. ,  572. ,  582. ,  592. ,  602. ,  612. ,
         622. ,  632. ,  642. ,  652. ,  662. ,  672. ,  682. ,  692. ,
         702. ,  712. ,  722. ,  732. ,  742. ,  752. ,  762. ,  772. ,
         782. ,  792. ,  802. ,  812. ,  822. ,  832. ,  842. ,  852. ,
         862. ,  872. ,  882. ,  892. ,  902. ,  912. ,  922. ,  932. ,
         942. ,  952. ,  962. ,  972. ,  982. ,  992. , 1002. , 1012. ,
        1022. , 1032. , 1042. , 1052. , 1062. , 1072. , 1082. , 1092. ,
        1102. , 1112. , 1122. , 1132. , 1142. , 1152. , 1162. , 1172. ,
        1182. , 1192. , 1202. , 1212. , 1217. , 1222. , 1227. , 1232. ,
        1237. , 1242. , 1247. , 1252. , 1257. , 1262. , 1267. , 1272. ,
        1277. , 1282. , 1287. , 1292. , 1297. , 1302. , 1307. , 1312. ,
        1317. , 1322. , 1327. , 1332. , 1337. , 1342. , 1347. , 1352. ,
        1357. , 1362. , 1367. , 1372. , 1377. , 1382. , 1387. , 1392. ,
        1397. , 1402. , 1407. , 1412. , 1417. , 1422. , 1427. , 1432. ,
        1437. , 1442. , 1447. , 1452. , 1457. , 1462. , 1467. , 1472. ,
        1477. , 1482. , 1487. , 1492. , 1497. , 1502. , 1507. , 1512. ,
        1517. , 1522. , 1532. , 1542. , 1552. , 1562. , 1572. , 1582. ,
        1592. , 1602. , 1612. , 1622. ])
    
    # min grid size
    min_grd_size = 0.05
    
    return LGR_sizes_xy, LGR_depths, min_grd_size

@pytest.fixture(scope='function')
def bbox_fixture():
    """ fixture for boundingbox 
    """

    bbox_casing = [0.762, 6, 21, 6, 21, 10, 16, 1, 16]
    bbox_cement_bond = [0.762, 6, 21, 6, 21, 10, 16]
    bbox_barrier = [0.3397, 10, 16, 10, 16, 14, 35]

    return bbox_casing, bbox_cement_bond, bbox_barrier
