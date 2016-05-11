input_file = "A1.csv"
MAXERROR = 10

"""
Script for ion selection for Joachm
Feb 22 2016
Author: Sean O'Callaghan
spoc@unimelb.edu.au
"""

from Tkinter import Tk, BOTH, Listbox, StringVar, END, Scrollbar, RIGHT, \
    VERTICAL, Y, TOP, LEFT, IntVar
from ttk import Frame, Label, Style, Entry, Button, Checkbutton


from numpy import arange, sin, pi
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt

# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure


from IonPickerEngine import deisotope, check_for_gaps


class Example(Frame):
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
         
        self.parent = parent        
        
      
        self.parent.title("Jo Ion Picker") 
        
        self.pack(fill=BOTH, expand=1)

        self.compounds = deisotope(filename=input_file)
        
        self.error_label = Label(self, text="Error in ppm")     
        self.error_label.place(x=20, y=730)

        self.entry = Entry(self, text='error in ppm')
        self.entry.place(x=20, y=750)
        self.entry.insert(10,'5')
        self.b = Button(self, text="ReCalc Error", width=15, \
                        command=self.callback)
        self.b.place(x=20, y=770)

        self.b_output = Button(self, text="Export", width=10, command=self.write_output)
        self.b_output.place(x=20, y=800)

        #self.b_output = Button(self, text="Allowed", width=10, command=self.only_allowed_ions)
        #self.b_output.place(x=20, y=830)


        self.gaps=IntVar()
        self.check_gaps = Checkbutton(self, text='Remove Gaps', variable=self.gaps,onvalue=1, offvalue=0, command=self.remove_gaps_call)
        #self.check.pack()
        self.check_gaps.place(x=20, y=830)


        self.scrollbar = Scrollbar(self, orient=VERTICAL)
        self.lb = Listbox(self, height=46, yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.lb.yview)
        self.scrollbar.pack(side=LEFT, fill=Y)

        
        for compound in self.compounds:
            print "found", compound.get_base_peak()
            mzs = compound.get_mz_list()
            num_mzs = len(mzs)
            entry = str(mzs[0]) + "    " + str(num_mzs)
            self.lb.insert(END, entry)
            
            
        self.lb.bind("<<ListboxSelect>>", self.onSelect)    
            
        self.lb.place(x=20, y=20)
        

        self.var = StringVar()
        #self.label = Label(self, text=0, textvariable=self.var)        
        #self.label.place(x=20, y=710)

        self.mz_label = Label(self, text="M/Z        Num Ions")
        self.mz_label.place(x=20, y=0)

        f = Figure(figsize=(8,11), dpi=100)
        self.ax = f.add_subplot(111)
        mass_list = self.compounds[0].get_mz_list()
        print mass_list
        intensity_list = self.compounds[0].get_intensity_list()
        mass_spec_plot = self.ax.bar(mass_list, intensity_list,\
        width=0.05)
        min_mz = mass_list[0]
        max_mz = mass_list[-1]
        self.ax.set_xlim([min_mz-1, max_mz+1])
        self.ax.set_ylim([0, 1.1*max(intensity_list)])
        self.ax.set_xticks(mass_list)
        self.ax.set_xticklabels(mass_list, rotation=45)
        self.ax.set_title("Base Ion:" + str(mass_list[0]))

        self.canvas = FigureCanvasTkAgg(f, master=self)
        self.canvas.show()
        self.canvas.get_tk_widget().pack(side=RIGHT)

    def onSelect(self, val):
      
        sender = val.widget
        idx = sender.curselection()
        value = sender.get(idx)   

        self.var.set(value)

        mz_to_search = value.split()[0]

        self.ax.clear()
        for i, compound in enumerate(self.compounds):
            if float(mz_to_search) == compound.get_base_peak():
                index = i
        mass_list = self.compounds[index].get_mz_list()
        print mass_list
        intensity_list = self.compounds[index].get_intensity_list()
        mass_spec_plot = self.ax.bar(mass_list, intensity_list,\
        width=0.05)
        min_mz = mass_list[0]
        max_mz = mass_list[-1]
        self.ax.set_xlim([min_mz-1, max_mz+1])
        self.ax.set_ylim([0, 1.1*max(intensity_list)])
        self.ax.set_xticks(mass_list)
        self.ax.set_xticklabels([str(mass) for mass in mass_list], rotation=45)
        self.ax.set_title("Base Ion:" + str(mass_list[0]))
        
        self.canvas.draw()

    def only_allowed_ions(self):

        ion_list = []

        fp = open('ions.csv', 'r')
        lines = fp.readlines()
        for line in lines:
            ion_list.append(float(line))
        #print ion_list
            
        self.compounds = deisotope(filename=input_file,max_error = float(self.entry.get()))

        new_compound_list = []
        

        for compound in self.compounds:
            for ion in ion_list:
                error_Da = float(self.entry.get())*compound.get_base_peak()/1e6
                if compound.get_base_peak() > ion-error_Da \
                   and compound.get_base_peak() < ion + error_Da:
                    new_compound_list.append(compound)

        for compound in new_compound_list:
            print "compound:",compound.get_base_peak()
        
        self.lb.delete(0, END)
        for compound in new_compound_list:
            mzs = compound.get_mz_list()
            num_mzs = len(mzs)
            entry = str(mzs[0]) + "    " + str(num_mzs)
            self.lb.insert(END, entry)
        
        

    def callback(self):
        self.compounds = deisotope(filename=input_file, max_error = float(self.entry.get()))
        self.lb.delete(0, END)
        for compound in self.compounds:
            mzs = compound.get_mz_list()
            num_mzs = len(mzs)
            entry = str(mzs[0]) + "    " + str(num_mzs)
            self.lb.insert(END, entry)
        print self.entry.get()


    def remove_gaps_call(self):
        self.compounds = deisotope(filename=input_file,max_error = float(self.entry.get()))
        check_for_gaps(self.compounds)
        self.lb.delete(0, END)
        for compound in self.compounds:
            mzs = compound.get_mz_list()
            num_mzs = len(mzs)
            entry = str(mzs[0]) + "    " + str(num_mzs)
            self.lb.insert(END, entry)
        print self.entry.get()


    def do_list_update(self):
        l_multi = self.multi.get()
        l_missing = self.missing.get()
        l_m1gtm2 = self.m1gtm2.get()
        #  possible situations: 000, 001,010,100, 011,101, 110, 111
        if l_multi == 0 and l_missing == 0 and l_m1gtm2 ==0:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                num_mzs = len(mzs)
                entry = str(mzs[0]) + "    " + str(num_mzs)
                self.lb.insert(END, entry)
                
        elif l_multi == 0 and l_missing == 0 and l_m1gtm2 ==1:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                num_mzs = len(mzs)
                intensities = compound.get_intensity_list()
                if intensities[-1] <= intensities[0]:
                    entry = str(mzs[0]) + "    " + str(num_mzs)
                    self.lb.insert(END, entry)
                    
        elif l_multi == 0 and l_missing == 1 and l_m1gtm2 ==0:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                num_mzs = len(mzs)
                if mzs[-1] - mzs[0] <1.75: # margin of error allowed here
                    entry = str(mzs[0]) + "    " + str(num_mzs)
                    self.lb.insert(END, entry)
                    
        elif l_multi == 1 and l_missing == 0 and l_m1gtm2 ==0:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                num_mzs = len(mzs)
                if num_mzs >1:
                    entry = str(mzs[0]) + "    " + str(num_mzs)
                    self.lb.insert(END, entry)
                    
        elif l_multi == 0 and l_missing == 1 and l_m1gtm2 ==1:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                intensities = compound.get_intensity_list()
                num_mzs = len(mzs)
                if mzs[-1] - mzs[0] <1.75 and intensities[-1] <= intensities[0]:
                    entry = str(mzs[0]) + "    " + str(num_mzs)
                    self.lb.insert(END, entry)
                
        elif l_multi == 1 and l_missing == 0 and l_m1gtm2 ==1:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                intensities = compound.get_intensity_list()
                num_mzs = len(mzs)
                if num_mzs >1 and intensities[-1] <= intensities[0]:
                    entry = str(mzs[0]) + "    " + str(num_mzs)
                    self.lb.insert(END, entry)
                
        elif l_multi == 1 and l_missing == 1 and l_m1gtm2 ==0:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                num_mzs = len(mzs)
                if num_mzs >1 and mzs[-1] - mzs[0] <1.75:
                    entry = str(mzs[0]) + "    " + str(num_mzs)
                    self.lb.insert(END, entry)
                
        elif l_multi == 1 and l_missing == 1 and l_m1gtm2 ==1:
            self.lb.delete(0, END)
            for compound in self.compounds:
                mzs = compound.get_mz_list()
                intensities = compound.get_intensity_list()
                num_mzs = len(mzs)
                if num_mzs >1 and mzs[-1] - mzs[0] <1.75 and \
                    intensities[1] <= intensities[0]:
                    entry = str(mzs[0]) + "    " + str(num_mzs)
                    self.lb.insert(END, entry)
        else:
            pass # error!
            

    def write_output(self):
        op = open('edited_output.csv', 'w')

        op.write('mz, intensity, mz, intensity, mz, intensity\n')
        items = self.lb.get(0,END)
        output_list = []
        for item in items:
            mz_val = item.split(' ')[0]
            for compound in self.compounds:
                if float(mz_val) == compound.get_base_peak():
                    mzs = compound.get_mz_list()
                    intensities = compound.get_intensity_list()
                    for i, mz in enumerate(mzs):
                        op.write(str(mz) + ',' + str(intensities[i]) + ',')
                    op.write('\n')
        op.close()
        

def main():
  
    root = Tk()
    ex = Example(root)
    root.geometry("1050x900+300+300")
    root.mainloop()  


if __name__ == '__main__':
    main()  
