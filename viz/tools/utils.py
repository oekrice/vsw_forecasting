"""
    Utilities.

    ary -- Oct 2022
"""
import os
import datetime
import time

# ------------------------------------------------------------------------------
def compile_f90(codepath="./viz/fortran/"):
    """
    Compile the fortran tracing code.
    """
    oldd = os.getcwd()
    os.chdir(codepath)
    os.system("make clean")
    os.system("make")
    os.chdir(oldd)


# ------------------------------------------------------------------------------
def list_snaps(path, prefix="b"):
    """
    List all netcdf snapshot files in given directory.
    """

    snaps = []
    for file in os.listdir(path):
        if file.startswith(prefix) and len(file) == len(prefix) + 14:
            snaps.append(file)
    snaps.sort()

    return snaps


# ------------------------------------------------------------------------------
def make_im_dir(path, name="im"):
    """
    Make a subdirectory in a run directory.
    """
    if not os.path.exists(path + name):
        os.makedirs(path + name)

#------------------------------------------------------------------------------
def toYearFraction(date):
    """
    Convert datetime object to fractional year.
    """
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime.datetime(year=year, month=1, day=1)
    startOfNextYear = datetime.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed / yearDuration

    return date.year + fraction
