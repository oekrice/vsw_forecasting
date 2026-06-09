import numpy as np
import sys
import os
import datetime


class Parameters(object):
    """
    Writer for Dumfric parameter (.nml) files.
    
    Attributes:
    ----------
      	     
      outputpath: String
      Path of output directory 
      
      starttime: String
      Start of simulation in format 'yyyymmdd.hh'
      
      endtime: String 
      End of simulation in format 'yyyymmdd.hh'
  
      nr: int
      Grid size in rho
      
      ns: int
      Grid size in s
      
      nph: int
      Grid size in phi (must be even for polar boundary condition)
      
      rss: float
      Outer boundary radius [in Rsun]
      
      cflfact: float
      CFL factor
      
      snapshot_cadence: int
            
      diagnostic_cadence: int
                        
      nu0: float
            
      ohmic_ratio: float
            
      ohmic_enhanced: float
            
      hyper_ratio: float     
      
      snapshot_variables: String
      
      driver_v: String
      '.TRUE.' or '.FALSE.'
      
      driver_em: String
      '.TRUE.' or '.FALSE.'
      
      cadence_em: int
      Duration of region emergence
    
  """

    def __init__(
        self,
        outputpath,
        starttime,
        endtime,
        nr=15,
        ns=30,
        nph=90,
        rss=2.5,
        cflfact=0.1,
        snapshot_cadence=24,
        diagnostic_cadence=1,
        nu0=0.36e-5,
        ohmic_ratio=0.0,
        ohmic_enhanced=0.2,
        hyper_ratio=1.2e-7,
        snapshot_variables="JBE",
        driver_v = '.TRUE.',
        driver_em = '.FALSE.',
        cadence_em = 24
    ):
        self.outputpath = outputpath
        self.nr = nr
        self.ns = ns
        self.nph = nph
        self.rss = rss
        self.starttime = starttime
        self.endtime = endtime
        self.outputpath = outputpath + "/"
        self.cflfact = cflfact
        self.snapshot_cadence = snapshot_cadence
        self.diagnostic_cadence = diagnostic_cadence
        self.nu0 = nu0
        self.ohmic_ratio = ohmic_ratio
        self.ohmic_enhanced = ohmic_enhanced
        self.hyper_ratio = hyper_ratio
        self.snapshot_variables = snapshot_variables
        self.driver_v = driver_v
        self.driver_em = driver_em
        self.cadence_em = cadence_em

        # Create output directory:
        print("Create output directory") 
        os.system("rm -r " + outputpath.replace(" ", "\ "))
        os.system("mkdir -p " + outputpath.replace(" ", "\ "))

    def write_prepared_file(self):
        """
      Write prepared.nml file for fortran

    """
        print("Writing prepared.nml")

        date = datetime.date.today().strftime("%d %B %Y")
        fid = open(self.outputpath + "prepared.nml", "w")
        fid.write("&prepared\n")
        fid.write("!\n")
        fid.write(
            "! DuMFric - auto-generated parameter file - created " + date + "\n"
        )
        fid.write("! -- DO NOT EDIT --\n")
        fid.write("!\n")
        fid.write("   start_time = " + self.starttime + ",\n")
        fid.write("   end_time = " + self.endtime + ",\n")
        fid.write("   nr = %i,\n" % self.nr)
        fid.write("   ns = %i,\n" % self.ns)
        fid.write("   np = %i,\n" % self.nph)
        fid.write("   rss = %g,\n" % self.rss)
        fid.write("   driver_v = " + self.driver_v + ",\n")
        fid.write("   driver_em = " + self.driver_em + ",\n")
        fid.write('   cadence_em = %i,\n' % self.cadence_em)
        fid.write("   /\n")
        fid.close()

    def write_user_file(self, filename="user.nml"):
        """
      Write user.nml file for fortran

    """
        print("Writing user.nml")

        fid = open(filename, "w")
        fid.write("&user\n")
        fid.write('datadir = "' + self.outputpath[:-1] + '",\n')
        fid.write("cflfact = %g,\n" % self.cflfact)
        fid.write("snapshot_cadence = %i,\n" % self.snapshot_cadence)
        fid.write("diagnostic_cadence = %i,\n" % self.diagnostic_cadence)
        fid.write('snapshot_variables = "' + self.snapshot_variables + '",\n')
        fid.write("nu0 = %g,\n" % self.nu0)
        fid.write("ohmic_ratio = %g,\n" % self.ohmic_ratio)
        fid.write("ohmic_enhanced = %g,\n" % self.ohmic_enhanced)
        fid.write("hyper_ratio = %g,\n" % self.hyper_ratio)
        fid.write("/\n")
        fid.close()
