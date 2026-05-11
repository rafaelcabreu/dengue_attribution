import sys

ok = 0
fail = 0

def check(label, fn):
    global ok, fail
    try:
        fn()
        print(f"  [PASS] {label}")
        ok += 1
    except Exception as e:
        print(f"  [FAIL] {label} — {e}")
        fail += 1

print("── Python package tests ─────────────────────────────")

# Import checks
check("import numpy",       lambda: __import__("numpy"))
check("import pandas",      lambda: __import__("pandas"))
check("import scipy",       lambda: __import__("scipy"))
check("import matplotlib",  lambda: __import__("matplotlib"))
check("import seaborn",     lambda: __import__("seaborn"))
check("import statsmodels", lambda: __import__("statsmodels"))
check("import geopandas",   lambda: __import__("geopandas"))
check("import shapely",     lambda: __import__("shapely"))
check("import unidecode",   lambda: __import__("unidecode"))

# Functional checks
def test_numpy():
    import numpy as np
    a = np.array([1.0, 2.0, 3.0])
    assert np.mean(a) == 2.0

def test_pandas():
    import pandas as pd
    df = pd.DataFrame({"x": [1, 2, 3], "y": [4, 5, 6]})
    assert df.shape == (3, 2)

def test_scipy():
    from scipy import stats
    result = stats.norm.cdf(0)
    assert abs(result - 0.5) < 1e-10

def test_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [4, 5, 6])
    plt.close(fig)

def test_statsmodels():
    import numpy as np
    import statsmodels.api as sm
    x = sm.add_constant(np.array([1.0, 2.0, 3.0, 4.0]))
    y = np.array([2.0, 4.0, 5.0, 4.0])
    sm.OLS(y, x).fit()

def test_shapely():
    from shapely.geometry import Point
    p = Point(0, 0).buffer(1)
    assert p.area > 3.1

def test_geopandas():
    import geopandas as gpd
    from shapely.geometry import Point
    gdf = gpd.GeoDataFrame({"val": [1]}, geometry=[Point(0, 0)], crs="EPSG:4326")
    assert len(gdf) == 1

def test_unidecode():
    from unidecode import unidecode
    assert unidecode("São Paulo") == "Sao Paulo"

check("numpy arithmetic",        test_numpy)
check("pandas DataFrame",        test_pandas)
check("scipy stats.norm.cdf",    test_scipy)
check("matplotlib plot (Agg)",   test_matplotlib)
check("statsmodels OLS",         test_statsmodels)
check("shapely buffer",          test_shapely)
check("geopandas GeoDataFrame",  test_geopandas)
check("unidecode São Paulo",     test_unidecode)

print("─────────────────────────────────────────────────────")
print(f"Results: {ok} passed, {fail} failed")

if fail > 0:
    sys.exit(1)
