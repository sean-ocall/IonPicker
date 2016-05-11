MAXERROR = 10

import math, os
import matplotlib.pyplot as plt

class Compound(object):

    def __init__(self, mz, intensity):
        self.__base_peak = mz
        self.__base_intensity = intensity
        self.__mz_list = []
        self.__mz_list.append(mz)
        self.__intensity_list = []
        self.__intensity_list.append(intensity)
        self.__is_isotope = False # use to tag isotopic peaks, remove later
        self.__is_of_interest = False
    
    def add_isotope(self, mz, intensity):
        
        self.__mz_list.append(mz)
        self.__intensity_list.append(intensity)

    def get_intensity_list(self):
        return self.__intensity_list

    def get_mz_list(self):
        return self.__mz_list

    def get_base_peak(self):
        return self.__base_peak

    def get_base_intensity(self):
        return self.__base_intensity

    def get_length(self):
        return len(self.__mz_list)

    def set_is_isotope(self, state):
        self.__is_isotope = state
        #print "changed",self.get_base_peak(), "to", state

    def get_is_isotope(self):
        return self.__is_isotope

    def set_is_of_interest(self, state):
        self.__is_of_interest = state

    def get_is_of_interest(self):
        return self.__is_of_interest

    def remove_records(self, remove_list):
        first_removed = remove_list[0]
        self.__mz_list = self.__mz_list[0:first_removed]
        self.__intensity_list = self.__intensity_list[0:first_removed]


def plot_ms(mass_list, intensity_list, plot_title=" "):
        
    """ 
    @summary: Plots the mass spec given a list of masses and intensities
        
    @param mass_list: List of the m/z values to be plotted
    @type mass_list: listType of floatType

    @param mass_list: List of the intensity values to be plotted
    @type mass_list: listType of floatType
        
    @param plot_title: A label for the plot
    @type plot_title: String Type
    
    @author: Sean O'Callaghan
    """
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
        
    # to set x axis range find minimum and maximum m/z channels
    max_mz = mass_list[-1]
    min_mz = mass_list[0]
        
    

    mass_spec_plot = plt.bar(mass_list, intensity_list,\
        width=0.05)
        
    x_axis_range = plt.xlim([min_mz-1, max_mz+1])
    y_axis_range = plt.ylim([0, 1.1*max(intensity_list)])
    x_axis_ticks = plt.xticks(mass_list)
    t = ax.set_title(plot_title)
    #plt.show()
    plt.savefig(plot_title +'.png', bbox_inches='tight')
    plt.close(fig)

#def filter_compounds(compound_list):

#    for compound in compound_list:
#        intensity_list = compound.get_intensity_list()
        

def deisotope(filename='jo-input.csv', output_file='output.csv',max_error=MAXERROR, check_gaps=True):
    # want to express error in mDa
    
    # we check how far away m/z values are from each other

    mass_C13 = 1.003355

    # then we divide by 1.034
    # We allow 3 isotopes at 100m/z, 6 at 2000m/z
    # Equation of line describing possible number of isotopes 
    # as mz increases
    # num_isotopes = (0.001578947)*mz + 2.842105263



    fp = open(filename, 'r')
    mzlist = []

    compound_list = []

    ion_list = []

    files = os.listdir(".")
    if 'ions.csv' not in files:
        print "*************************************************"
        print "*                                               *"
        print "*         Error: ions.csv not present!!!        *"
        print "*                                               *"
        print "*************************************************"

    ions_fp = open('ions.csv', 'r')
    lines = ions_fp.readlines()
    for line in lines:
        ion_list.append(float(line))
        print float(line)
    
    for mz_intensity in fp:
        #print mz
        mz, intensity = mz_intensity.split(',')
        if mz_intensity[0:3] == 'm/z':
            continue
        error_Da = max_error*float(mz)/1e6
        for ion in ion_list:
            #print ion-error_Da, ":",mz,":",ion+error_Da
            if float(mz) > ion-error_Da and float(mz) <= ion+error_Da:
                compound = Compound(float(mz), float(intensity.strip('\n')))
                #print "Found",compound.get_base_peak()
                compound.set_is_of_interest(True)
                #print compound.get_base_peak(), "set to True"
                compound_list.append(compound)
            elif float(mz) > ion+error_Da and float(mz) < ion+error_Da+10:
                compound = Compound(float(mz), float(intensity.strip('\n')))
                compound.set_is_isotope(True)
                #print "found", compound.get_base_peak(),
                #if compound.get_is_isotope():
                #    print "Is Isotope (maybe)?"
                
                #print "created compound:", float(mz), float(intensity.strip('\n'))
                compound_list.append(compound)

    print "Compound list length:", len(compound_list)

    for i,compound in enumerate(compound_list):
        if not compound.get_is_isotope():
            
            mz = compound.get_base_peak()
            for j, later_compound in enumerate(compound_list[i+1:]):
                later_mz = later_compound.get_base_peak()
                max_isotopes = 10#(0.001578947)*mz + 2.842105263 # see top for details

                #Convert error from ppm to Dalton
                error_Da = max_error*mz/1e6

                # if it's not very close
                # and if it's less than max_error away from an integer multiple of
                # 1.034 (C13 mass)
                # and it's less than max_isotopes away
                if (later_mz - mz) > 1.0 and \
                        math.fabs(mass_C13*(round(later_mz-mz))-(later_mz-mz))\
                        < error_Da \
                        and mass_C13*(round(later_mz-mz))<max_isotopes:
                    #print "%f,%d,%f"%(later_mz, j+1, mz)
                    compound.add_isotope(later_mz, later_compound.get_base_intensity())
                    #later_compound.set_is_isotope(True)

    lines = []
   
    if check_gaps == True:
        check_for_gaps(compound_list)
   
    # See how many lines we need
    line_counter = 0
    for compound in compound_list:
        if not compound.get_is_isotope():
            if len(compound.get_mz_list()) > line_counter:
                line_counter = len(compound.get_mz_list())

    op = open(output_file,'w')
    
    

    # Create blank lines
    for i in range(line_counter):
        lines.append("")
        
    line_1 = ''
    for compound in compound_list:
        if not compound.get_is_isotope():
            line_1 = line_1+'m/z, intensity,'
            mz_list = compound.get_mz_list()
            intensity_list = compound.get_intensity_list()

            for i, line in enumerate(lines):
                if i < len(mz_list):
                    newline = line + str(mz_list[i]) + ',' + \
                           str(intensity_list[i]) + ','
                    lines[i] = newline
                else:
                    newline = line + ',,'
                    lines[i] = newline

    op.write(line_1 + '\n')
    for line in lines:
        op.write(line + '\n')


    fp.close()
    op.close()

    mzs = []

    for compound in compound_list:
        if not compound.get_is_isotope():
            mzs.append(str(compound.get_base_peak()))

    new_compound_list = []
    for compound in compound_list:
        if compound.get_is_of_interest():
            #print compound.get_base_peak(), "sent to GUI"
            new_compound_list.append(compound)

    #check_for_gaps(new_compound_list)

    return new_compound_list


def check_for_gaps(compounds):

    for compound in compounds:
        remove_list = []
        mz_list = compound.get_mz_list()
        for i, mz in enumerate(compound.get_mz_list()[:-1]):
            if mz_list[i+1] - mz >1.2:
                remove_list.append(i+1)

        if len(remove_list) >0:
            compound.remove_records(remove_list)
        
    



