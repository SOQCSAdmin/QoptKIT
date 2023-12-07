###############################################################################
# QOPTKIT interface with Python. Python side.                                   #  
#                                                                             #
#                                                                             #
# Copyright 2023 Maynooth University                                          #
#                                                                             #
# Licensed under the Apache License, Version 2.0 (the "License");             #
# you may not use this file except in compliance with the License.            #
# You may obtain a copy of the License at                                     #
#                                                                             #
#     http://www.apache.org/licenses/LICENSE-2.0                              #
#                                                                             #
# Unless required by applicable law or agreed to in writing, software         #
# distributed under the License is distributed on an "AS IS" BASIS,           #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
# See the License for the specific language governing permissions and         #
# limitations under the License.                                              #
###############################################################################

import sys
import os
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from ctypes import cdll,c_char,c_int,c_long,c_double,POINTER
from qiskit import QuantumCircuit, transpile  
from qiskit.extensions import UnitaryGate


#------------------------------------------------------------------------------#      
# CPP library configuration
#------------------------------------------------------------------------------#      
cpplib='libqoptkit.so'

#------------------------------------------------------------------------------#      
# Open QOPTKIT C++ interface                                                   #
#------------------------------------------------------------------------------#      
qoptkit = cdll.LoadLibrary(os.path.dirname(__file__) + '/' + cpplib)


#------------------------------------------------------------------------------#     
# Memory management auxiliary methods.                                         #
# Free the memory of an object created in C++                                  #
# Not valid for classes                                                        #
#------------------------------------------------------------------------------#      
def free_ptr(array):
    """

    Free memory assigned to pointers defined in C++
    
    :array(POINTER): Vector pointer
    
    """    
    array.restype=POINTER(c_char)
    qoptkit.free_ptr(array)

#------------------------------------------------------------------------------#      
# Converter MTX->int*                                                          #
#------------------------------------------------------------------------------#      
def to_int_ptr(mtx):
    """

    Converts a mtx into a integer pointer.
    
    :mtx(list[][]): Matrix
    
    """  
    n=len(mtx)
    m=len(mtx[0])
    length=n*m


    send=(c_int*length)()
    for i in range(0,n):
        for j in range(0,m):            
            send[i*m+j]=mtx[i][j]

    return send,n,m


#------------------------------------------------------------------------------#
# Wrapper for the C++ QOPTKIT configuration method                                   #
# Configuration of the default maximum number of photons of the simulation     #
#------------------------------------------------------------------------------#     
def cfg_qoptkit(nph):
    """

    It defines the default maximum number of photons to be used by qoptkit.
    
    :nph (int): Maximum number of photons used by qoptkit.
    
    """  
    qoptkit.all_cfg_qoptkit(nph);   



#------------------------------------------------------------------------------#      
# Class gcircuit: Graphical circuit                                            #
# Auxiliary class implementing the graphical representation of                 #
# a Quantum Optical Device.                                                    #
#------------------------------------------------------------------------------#      
class gcircuit(object):
    #---------------------------------------------------------------------------      
    # Creates a gcircuit object
    #---------------------------------------------------------------------------      
    def __init__(self):              
        self.list=[]    # List
        self.nl=0       # Numbe of elements in the list
        self.newrow=0   # Number of calls to newline. Affects the plot dimensions.
      
    #---------------------------------------------------------------------------      
    # Adds a gate to the plot specifying input channels.
    # The gate is plotted with a custom text.
    #---------------------------------------------------------------------------            
    def add_gate(self, chlist, text):
        nch=len(chlist)

        
        # Rewire
        ch1=chlist[0]
        for i in range(0,nch):
            ch2=chlist[i]
            if(ch2!=ch1+i):
                self.rewire(ch2,ch1+i)    
            
        # Plot gate
        self.element(ch1,nch,text, False)
        
        # Rewire
        for i in range(0,nch):
            ch2=chlist[i]
            if(ch2!=ch1+i):
                self.rewire(ch1+i,ch2)  
                        
    #---------------------------------------------------------------------------      
    # Defines a wire in a channel with length depth.
    #---------------------------------------------------------------------------      
    def add_wire(self,ch,depth):   
        self.list.append(['add_wire',ch,depth])
        self.nl=self.nl+1
        
    #---------------------------------------------------------------------------      
    # Defines a circuit element
    #---------------------------------------------------------------------------      
    def element(self,ch,nch,text, colored=False): 
        self.list.append(['element',ch,nch,text, colored])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Defines a rewire between two channels
    #---------------------------------------------------------------------------      
    def rewire(self,ch1,ch2,initial=False): 
        self.list.append(['rewire',ch1,ch2,initial])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Defines a single channel circuit element
    #---------------------------------------------------------------------------      
    def one_gate(self,ch,text, colored=False): 
        self.element(ch,1,text, colored)        

    #---------------------------------------------------------------------------              
    # Defines a two channel circuit elements
    # Rewires for non consecutive channels are calculated automatically
    #---------------------------------------------------------------------------      
    def two_gate(self,ch1,ch2,text, colored=False): 
        if(ch2!=ch1+1):
            self.rewire(ch2,ch1+1)
        
        self.element(ch1,2,text, colored)
    
        if(ch2!=ch1+1):
            self.rewire(ch1+1,ch2)

    #---------------------------------------------------------------------------      
    # Defines a separator
    #---------------------------------------------------------------------------      
    def separator(self): 
        self.list.append(['separator'])
        self.nl=self.nl+1

    #---------------------------------------------------------------------------      
    # Defines a newline
    #---------------------------------------------------------------------------      
    def newline(self): 
        self.list.append(['newline'])
        self.nl=self.nl+1
        self.newrow=self.newrow+1
        
    #---------------------------------------------------------------------------              
    # Draws an empty wire in a channel for a given depth
    #---------------------------------------------------------------------------      
    def g_add_wire(self,ch,depth):
        if(self.ended[ch]==0):
            init=self.display[ch];
            end=depth
        
            if(depth>init):
                height=0.5+(1+self.pad)*ch+self.line*self.width_line            
                if init>0:                               
                    plt.plot([2*init, 2*end],[height, height], color='black', linewidth=4)
                else:
                    plt.plot([0.25,2*end],[height, height], color='black', linewidth=4)
                    text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
                    plt.text(0.0,height,str(ch),**text_kwargs)
                
                self.display[ch]=depth

    #---------------------------------------------------------------------------      
    # Draws a representation of a general circuit element
    #---------------------------------------------------------------------------      
    def g_element(self,ch,nch,text, colored=False): 
        maxdepth=0
        for i in range(ch,ch+nch):
            maxdepth=max(maxdepth,self.display[i])
            self.ended[i]=0
            
        init=maxdepth
        end=init+1;    
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch];
            end=init+1;
        
        for i in range(ch,ch+nch):
            self.g_add_wire(i,init)
            self.display[i]=end                                
            
        height=0.5+(1+self.pad)*ch+self.line*self.width_line                        
        add=(nch-1)      

        color='mediumpurple'
        if colored==True:
            color='mediumorchid'            
        rect = patches.Rectangle(((2*init+0.2), (height-0.5+0.2)),
                                  (2*end-0.2-(2*init+0.2)), 
                                  ((height+add+0.5-0.2)-(height-0.5+0.2)),
                                 linewidth=4, edgecolor='black', facecolor=color)
        self.ax.add_patch(rect)

        
        for i in range(0,nch):
            heightch=0.5+(1+self.pad)*(ch+i)+self.line*self.width_line            
            
            plt.plot([2*init,(2*init+0.2)],[heightch, heightch], color='black', linewidth=4)
            plt.plot([(2*end-0.2),2*end],[heightch, heightch], color='black', linewidth=4)
                        
        
        text_kwargs = dict(ha='center', va='center',multialignment='left', fontsize=self.size_font-2, color='black')
        plt.text((2*init+1),(height+(nch-1)/2),str(text),**text_kwargs)

    #---------------------------------------------------------------------------      
    # Draws a rewire between two channels
    #---------------------------------------------------------------------------      
    def g_rewire(self,ch1,ch2,initial): 
        maxdepth=0
        chi=min(ch1,ch2)
        cho=max(ch1,ch2)
        for chaux in range(chi,cho+1):
            maxdepth=max(maxdepth,self.display[chaux])
            
        if maxdepth>0:
            init=maxdepth
            end=init+1;
            if end>(self.depth-2):
                self.g_newline()
                init=self.display[ch1];
                end=init+1;

            self.g_add_wire(chi,init)   
            self.display[chi]=end        
            self.g_add_wire(cho,init)   
            self.display[cho]=end        
            for chaux in range(chi+1,cho):
                self.g_add_wire(chaux,end)   
                self.display[chaux]=end
                                
            height1=0.5+(1+self.pad)*chi+self.line*self.width_line                       
            height2=0.5+(1+self.pad)*cho+self.line*self.width_line                       
        
            if self.ended[chi]==0:
                plt.plot([2*init,2*end],[height1, height2], color='black', linewidth=4)
            if self.ended[cho]==0:
                plt.plot([2*init,2*end],[height2, height1], color='black', linewidth=4)   
        
        aux=self.ended[chi]
        self.ended[chi]=self.ended[cho]
        self.ended[cho]=aux

    #---------------------------------------------------------------------------      
    # Draws the initial channel lines
    #---------------------------------------------------------------------------  
    def g_init_ph(self,ch,n):    
        self.ended[ch]=0
        self.nph[ch]=self.nph[ch]+n
        height=0.5+(1+self.pad)*ch   
        
        init=self.display[ch]
        end=init;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch]
            end=init
        
        self.display[ch]=end 
        posx=2*init
         
                        
        plt.plot([posx+0.25,posx+2],[height, height], color='black', linewidth=4)
        plt.plot([posx+1.5,posx+2],[height, height], color='black', linewidth=4)
        text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
        plt.text(posx+0.0,height,str(ch),**text_kwargs)   
        
        if(self.display[ch]<=0):
            self.display[ch]=1
            
    #---------------------------------------------------------------------------      
    # Draws the ending channel lines
    #---------------------------------------------------------------------------     
    def g_final_dec(self,ch,post,pol):        
        self.ended[ch]=1
        init=self.display[ch]
        end=init+1;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[ch]
            end=init+1



        self.display[ch]=end
        height=0.5+(1+self.pad)*ch+self.line*self.width_line                        
        
        text_kwargs = dict(ha='center', va='center', fontsize=self.size_font, color='black')
        plt.text(2*init+0.25,height,str(ch),**text_kwargs)   
                   
    #---------------------------------------------------------------------------              
    # Draws a separator to clarify the plot
    #---------------------------------------------------------------------------      
    def g_separator(self, muted=0): 
        maxdepth=0
        for i in range(0,self.nch):
            maxdepth=max(maxdepth,self.display[i])

        init=maxdepth
        end=init+1;
        if end>(self.depth-2):
            self.g_newline()
            init=self.display[0];
            end=init+1;
               
        for i in range(0,self.nch):
            self.g_add_wire(i,end)
            self.display[i]=end
            
        if muted==0:
            plt.plot([2*init+1, 2*init+1],
                     [0+self.line*self.width_line , self.nch+self.line*self.width_line], 
                     color='gray',linestyle='dashed', linewidth=6)

    #---------------------------------------------------------------------------              
    # Draws a "hard" separator to identify changes of line
    #---------------------------------------------------------------------------      
    def g_hard_separator(self, width): 
        maxdepth=0
        for i in range(0,self.nch):
            maxdepth=max(maxdepth,self.display[i])

        init=maxdepth
        end=init+1;
               
        for i in range(0,self.nch):
            self.g_add_wire(i,end)
            self.display[i]=end
            
        plt.plot([2*init+1, 2*init+1],
                 [0+self.line*self.width_line , self.nch+self.line*self.width_line], 
                 color='black',linestyle='dashdot', linewidth=width)

    #---------------------------------------------------------------------------              
    # Advances the printing row.
    #---------------------------------------------------------------------------           
    def g_newline(self): 
        self.g_hard_separator(3)
        self.display=np.zeros(self.nch)
        self.line=self.line+1
        self.g_hard_separator(3)
        for i in range(0,self.nch):
            self.g_add_wire(i,1)
                    
    #---------------------------------------------------------------------------          
    # Plots the circuit
    #---------------------------------------------------------------------------      
    def show(self, nch, depth, sizexy, font, rows=1):
  
        dpi=sizexy;
        
        self.nch=nch
        self.size_font=int(font)        
        self.display=np.zeros(nch)
        self.nph=np.zeros(nch)
        self.pad=0.1
        self.line=0
        self.width_line=((sizexy+sizexy/10)*(nch-1)+sizexy)/dpi+0.5*sizexy/dpi
        self.depth=depth
        if rows>1:
            self.newrow=rows-1
        self.ended=[0]*nch
        
        self.fig, self.ax= plt.subplots()        
        self.fig.set_dpi(dpi)
        self.fig.set_size_inches( ((depth-1)*2*sizexy+sizexy)/dpi, (self.newrow+1)*((sizexy+sizexy/10)*nch+sizexy)/dpi, forward=True)
        self.ax.set_facecolor("white")
        plt.xlim([0, ((depth-1)*2*sizexy-0.8*sizexy)/dpi])
        plt.ylim([-0.5*sizexy/dpi, (self.newrow+1)*self.width_line])
        plt.gca().invert_yaxis()
        plt.axis('off')
        
        
        for i in range(0,nch):
            self.g_init_ph(i,-1)
            
        for i in range(0,self.nl):
            gate=self.list[i]
            
            if(gate[0]=='add_wire'):             
                self.g_add_wire(gate[1],gate[2])

            if(gate[0]=='element'): 
                self.g_element(gate[1],gate[2],gate[3], gate[4])

            if(gate[0]=='rewire'): 
                self.g_rewire(gate[1],gate[2],gate[3])
          

            if(gate[0]=='separator'): 
                self.g_separator(muted=0)
                
            if(gate[0]=='newline'): 
                self.g_newline()

        self.g_separator(muted=1)
        for i in range(0,nch):
            self.g_final_dec(i,-1,-1)
        plt.draw()
        plt.show()

    
#------------------------------------------------------------------------------#      
# Definition of a quantum optical circuit                                      #
#------------------------------------------------------------------------------#     
class qocircuit(object):
    """

    A quantum circuit consists in a set of interconnected quantum optical elements.

    :nch (int): Number of channels
    :dummy(automatic(bool)): If True it creates a dummy empty class. (Option used internally by the library).
    
    """

    #---------------------------------------------------------------------------      
    # Create circuit
    #---------------------------------------------------------------------------     
    def __init__(self, nch, dummy=False):
            
        func=qoptkit.qoc_new_qocircuit
        func.restype=c_long    
        if(dummy==False):
            self.obj = func(nch)
            self.disp=gcircuit()
            self.nch=nch

    #---------------------------------------------------------------------------              
    # Delete circuit        
    #---------------------------------------------------------------------------      
    def __del__(self):       
        qoptkit.qoc_destroy_qocircuit(c_long(self.obj))

    #---------------------------------------------------------------------------              
    # Return nunmber of channels          
    #---------------------------------------------------------------------------      
    def num_ch(self):
        """

        Returns the total number of channels of the circuit.

        :return(int): Number of channels of the circuit.
        """
        return qoptkit.qoc_num_levels(c_long(self.obj))


    #---------------------------------------------------------------------------      
    # Adds a beamsplitter to the circuit
    #---------------------------------------------------------------------------      
    def bs(self, ch1, ch2, theta, phi):
        """
            
        Adds a beamsplitter to the quantum circuit attached to channels ch1 and ch2.
        
        :ch1 (int): Beamsplitter input channel 1.
        :ch2 (int): Beamsplitter input channel 2.
        
        """
        # Logic
        qoptkit.qoc_beamsplitter(c_long(self.obj),c_int(ch1),c_int(ch2),c_double(theta),c_double(phi))
        
        # Graphics
        str1="     BS\n \u03B8="
        str2=str(round(theta,2))
        str3="º\n \u03C6="
        str4=str(round(phi,2))
        str5="º"
        strf=str1+str2+str3+str4+str5
        self.disp.two_gate(ch1, ch2,strf, False)

    #---------------------------------------------------------------------------              
    # Adds a phase shifter to the circuit
    #---------------------------------------------------------------------------      
    def ph(self, ch, phi):
        """
            
        Adds a phase shifter to the circuit in channel ch.
        
        :ch (int): Phase shifter input channel.
        :phi (float): Angle phi in degrees.
        
        """
        # Logic
        qoptkit.qoc_phase_shifter(c_long(self.obj),c_int(ch),c_double(phi))

        # Graphics
        str1="PS \u03C6="
        str2=str(round(float(phi),2))
        str3="º"
        strf=str1+str2+str3
        self.disp.one_gate(ch,strf,False)        


    #---------------------------------------------------------------------------      
    # Swaps to channels
    #---------------------------------------------------------------------------      
    def swp(self, ch1, ch2):
        """
            
        Adds a swap gate between two channels
        
        :ch1 (int): Channel 1.
        :ch2 (int): Channel 2.
        
        """
        # Logic
        qoptkit.qoc_rewire(c_long(self.obj),ch1,ch2)
        
        # Graphics
        self.disp.rewire(ch1, ch2) 

        
    #---------------------------------------------------------------------------      
    # Adds a gate using other circuit as the gate definition
    #---------------------------------------------------------------------------    
    def gate(self, chlist, qoc, text):
        """

        Adds a gate using other circuit as the gate definition.
   
        :chlist (list[int]): List of channels to which the new gate is attached
        :qoc (qocircuit): Circuit defining the gate
    
        """
        # Logic 
        length=len(chlist)
        send=(c_int*length)()
        for i in range(0,length):
            send[i]=chlist[i]

        qoptkit.qoc_add_gate(c_long(self.obj),send,c_long(qoc.obj))
        
        

        # Graphics
        self.disp.add_gate(chlist,text)

    #---------------------------------------------------------------------------      
    # Adds a graphical separator between two circuit sections
    #---------------------------------------------------------------------------    
    def barrier(self):
        """
            
        Adds a graphical separator between two circuit sections
        
        
        """        
        self.disp.separator()
       
    #---------------------------------------------------------------------------              
    # Plots the circuit
    #---------------------------------------------------------------------------      
    def show(self,depth=10,sizexy=100,font=18, rows=1):
        """
        
        Plots the quantum circuit.
        
        :depth (optional[int]): Number of plot layers. The longer the circuit the more needed.
        :sizexy(optional[int]): How many pixels by layer. Controls the size of the plot.
        :font (optional[int]):  Font size of the labels.
        :rows(int): Numbers of printing lines to be used by the plotter. |br| 
        
        """    
        self.disp.show(self.nch,depth,sizexy,font, rows)        

#------------------------------------------------------------------------------#                                              #
# Quantum photonic state                                                       #
#------------------------------------------------------------------------------#      
class state(object):
    """

    Quantum photonic state definition. A quantum state is a sum of kets multiplied by amplitudes.

    :level (int): Number of channels to describe the state.
    :st(optional(list[[complex,list[int]]]): List of kets described as a 2.tuple of a comlpex number representing the amplitude and a list of occupations by channel.
    :nph(optional(int)): Maximum number of photons of the state.
    :maxket(optional(int)): Maximum number of kets. (Internal memory).
    :dummy(automatic(bool)): If True it creates a dummy empty class. (Option used internally by the library).
    
    """   
    #---------------------------------------------------------------------------      
    # Create a state
    #---------------------------------------------------------------------------      
    def __init__(self, level,st=[],nph=-1,maxket=50,dummy=False):
                
        func=qoptkit.st_new_state
        func.restype=c_long
        if(dummy==False):
            self.obj = func(nph,level,maxket)

            nkets=len(st)
            if nkets>0:    
                for i in range(0,nkets):
                    ket=st[i]
                    ampl=ket[0]
                    occ=ket[1]
                    self.add_ket(ampl,occ)
            

    #---------------------------------------------------------------------------              
    # Delete a state
    #---------------------------------------------------------------------------      
    def __del__(self):
        qoptkit.st_destroy_state(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Calculate the braket coefficient between two states
    #---------------------------------------------------------------------------      
    def braket(self,state):   
        """

        Performs the bra-ket operation <bra|state> using the complex conjugate of this state as bra.
        
        :state(state): State in the right hand side of the braket operation.         
        :return (complex): The complex number result of the projection.
            
        """  
        func=qoptkit.st_braket
        func.restype=POINTER(c_double)
        array_ptr=2*c_double
        array_ptr=func(c_long(self.obj),c_long(state.obj))
        value=complex(c_double(array_ptr[0]).value,c_double(array_ptr[1]).value)        
        free_ptr(array_ptr)
        return value

    #---------------------------------------------------------------------------      
    # Normalizes a state
    #---------------------------------------------------------------------------      
    def normalize(self):
        """

        Normalizes the state.         
        
        """  
        qoptkit.st_normalize(c_long(self.obj))

        
    #---------------------------------------------------------------------------      
    #  Adds a new term (amplitude + ket) to a state directly from level definition
    #---------------------------------------------------------------------------      
    def add_ket(self, ampl, vec):   
        """

        Adds a new ket to the state using a vector with the level occupation.
         
        :ampl(complex):  Amplitude of the new ket.
        :vec(list[int]): List with the occupation of each channel in the new ket.
        
        """  
        length=len(vec)
        send=(c_int*length)()
        for i in range(0,length):
            send[i]=vec[i]
            
        qoptkit.st_add_term(c_long(self.obj),c_double(ampl.real), c_double(ampl.imag), send) 

    #---------------------------------------------------------------------------      
    # Post-selection by a projector. Auxiliary method. Inteded for internal use.
    #---------------------------------------------------------------------------      
    def prj_post_selection(self, prj):             
        func=qoptkit.st_post_selection
        func.restype=c_long
        obj=func(c_long(self.obj), c_long(prj.obj))
        newstate=state(1,dummy=True) 
        newstate.obj=obj
        return newstate
    
    #---------------------------------------------------------------------------      
    # Post-selection by a projector
    #---------------------------------------------------------------------------      
    def post_selection(self, cond):   
        """

        Post-selection over states by the condition defined by the list cond.
         
        :cond(list[int]): Description of the post-selection condition a list of occupations by channel. Channels with negative values of the occupation are not post-selected.
        :return (state): Post-selected state.
        
        """            
        
        prj=projector(self.num_ch(), cond)
        return self.prj_post_selection(prj)
    

    #---------------------------------------------------------------------------      
    # Encode a state into qubit encoding. (Path encoding)
    # Those kets that can not be encoded are ignored and the result is normalized
    # WARNING: It is responsibility of the library user to make sure that no ket
    # repetitions arise after encoding.
    #---------------------------------------------------------------------------      
    def encode(self, qmap,  qoc):
        """
        Encode a state into a qubit encoding. Those kets that can not be encoded are dismissed and the result is normalized. |br|
    
        **Warning!** t is responsibility of the library user to make sure that no ket repetitions arise after encoding.


        :qmap(list[][]): nx2 matrix with the qubit definitions. Each row has two entries specifying the channels that define the qubit.
        :qoc(qocircuit): Circuit to which the state is related.
        :return output(state): An encoded state.
        
        """  
        param=to_int_ptr(qmap)
        func=qoptkit.st_encode
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[1],c_long(qoc.obj))  
        aux=state(1,dummy=True)
        aux.obj=obj
        return aux 

   #----------------------------------------------------------------------------      
   # Decode a state from a qubit encoding into photon modes. (Path encoding)
   #----------------------------------------------------------------------------
    def decode(self, qmap, ancilla, qoc):
        """
        Decode a state from a qubit encoding into photon encoding. 
        Note that we need to define the values of the extra ancilla channels that are not included into the qubit encoding
        but they are needed to the circuit to work in the intended way


        :qmap(list[][]): nx2 matrix with the qubit definitions. Each row has two entries specifying the channels that define the qubit.
        :ancilla(list[int]): List with the values of the ancilla channels (from smaller to large channel number).
        :qoc(qocircuit): Circuit which the state is related.
        :return output(state): A decoded state.
        
        """  
        
        nq=len(qmap)
        nanz=len(ancilla)
        nch=2*nq+nanz   

      
        # Compute ancillas
        aux=[1]*nch
        for i in range(0,nq):
            aux[qmap[i][0]]=0
            aux[qmap[i][1]]=0
            
        k=0
        occ=[]
        for i in range(0,nch):
            if aux[i]==0:
                occ.append(0)
            if aux[i]==1:
                occ.append(ancilla[k])
                k=k+1
                
        anzstate= state(nch,maxket=1)
        anzstate.add_ket(1.0, occ)
        
        
        param=to_int_ptr(qmap)
        func=qoptkit.st_decode
        func.restype=c_long   
        obj=func(c_long(self.obj),param[0],param[1],c_long(anzstate.obj),c_long(qoc.obj))  
        aux=state(1,dummy=True)
        aux.obj=obj
        return aux 

    
    #---------------------------------------------------------------------------      
    # Print state
    #---------------------------------------------------------------------------      
    def prnt_state(self):
        """
            Prints a state.
        """  
                                    
        qoptkit.st_prnt_state(c_long(self.obj))
        sys.stdout.flush()
                                    

    #---------------------------------------------------------------------------          
    # Return the number of kets
    #---------------------------------------------------------------------------      
    def nkets(self):
        """
            Returns the number of kets stored.        
    
            :return (int): Number of kets.
        """
        return qoptkit.st_nkets(c_long(self.obj))
    
    #---------------------------------------------------------------------------          
    # Return the number of channels
    #---------------------------------------------------------------------------      
    def num_ch(self):
        """
            Returns the number of channels.
    
            :return (int): Number of channels.
        """
        return qoptkit.st_nlevels(c_long(self.obj))    

    #---------------------------------------------------------------------------      
    # Return the amplitude of probability of a ket
    #---------------------------------------------------------------------------      
    def campl(self, index):
        """
        
        Returns the amplitude  of probability of the required ket
    
        :index (int):  Index of the ket.
        :return (int): Amplitude of probability of the referred ket.
        
        """
        func=qoptkit.st_campl
        func.restype=POINTER(c_double)
        array_ptr=2*c_double
        array_ptr=func(c_long(self.obj),index)
        value=complex(c_double(array_ptr[0]).value,c_double(array_ptr[1]).value)        
        free_ptr(array_ptr)
        return value
          
    #---------------------------------------------------------------------------      
    # Return tag of the ket given its index.
    #---------------------------------------------------------------------------      
    def tag(self, index):
        """
        
        Returns the tag of the ket required.
    
        :index (int):  Index of a ket.
        :return (string): Occupation in string format.
        
        """
        nlevel=self.num_ch()
        func=qoptkit.st_tag
        func.restype=POINTER(c_char)
        array_ptr=c_char*nlevel
        array_ptr=func(c_long(self.obj),index)
        array=[0]*nlevel
        for i in range(0, nlevel):
            array[i]=array_ptr[i].decode('UTF-8')
        free_ptr(array_ptr)
        return "".join(array)
    

    #---------------------------------------------------------------------------      
    # Return state as a list of 2-uples of amplitudes and occupations.
    #---------------------------------------------------------------------------       
    def ket_list(self):
        """
        
        Return state as a list of 2-tuples of amplitudes and occupations.
    
        :return (list[[complex,[int]]): List of 2-tuples of a ket amplitude and a list of occupations by channel.
        
        """        
        nkets=self.nkets()
        st_list=list()
        for i in range(0,nkets):            
            st_list.append([self.campl(i),self.tag(i)])        
        
        return st_list
            
        
    
    #---------------------------------------------------------------------------      
    # Plots the probabilities of each ket graphically
    #---------------------------------------------------------------------------      
    def show(self, pmax=-1.0,sizex=8, sizey=5, dpi=100, angle=70, font=14):
        """
        
            Plots the outcomes in a bar diagram.
        
            :pmax(optional[float]): Maximum value of the probabiilty.
            :sizex(optional[float]): Size of the plot in the x direction in inches.
            :sizey(optional[float]): Size of the plot in the y direction in inches.
            :dpi(optional[float]): Density of points  by inch.
            :angle(optional[float]): Angle of the horizontal axis labels.
            :font(optional[fint]): Size of the font of the horizontal axis labels.
            
        """  
        nkets=self.nkets()
        taglist=list()
        problist=list()
        for i in range(0,nkets):
            taglist.append(self.tag(i))
            ampl=self.campl(i)
            problist.append(np.real(ampl*np.conjugate(ampl)))
  
        zipped_lists = zip(taglist, problist)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        staglist, sproblist = [ list(tuple) for tuple in  tuples]
        
        fig, ax = plt.subplots(1, 1)
        fig.set_size_inches(sizex,sizey)
        fig.set_dpi(dpi)
        plt.xticks(rotation = angle)
        if pmax>0.0:
            plt.ylim(0,pmax)
        
        y_pos = np.arange(len(staglist))
        pbar=ax.bar(y_pos, sproblist, align='center', alpha=0.5)
        ax.bar_label(pbar,fmt='%.3f')
        ax.tick_params(axis='x', which='major', labelsize=font)
        ax.set_axisbelow(True)
        plt.grid(visible=True, which='major', axis='y', color='0.8', linestyle='--', linewidth=2)
        plt.xticks(y_pos, staglist)
        plt.ylabel('P',fontsize=font)
        plt.title('Output')
       
        plt.draw()
        plt.show()
        
#---------------------------------------------------------------------------                  
# Wrapper for C++ QOPTKIT class projector
# Definition of a projector to perform post-selection. 
# Class intended for internal use.
#---------------------------------------------------------------------------      
class projector(object):
 
    #---------------------------------------------------------------------------      
    # Create a projector
    #---------------------------------------------------------------------------      
    def __init__(self, level,cond=[],nph=-1,maxket=50,dummy=False):        

                
        func=qoptkit.prj_new_projector
        func.restype=c_long
        if(dummy==False):
            self.obj = func(nph,level,maxket)

            
            nkets=len(cond)
            if nkets>0:    
                self.set_cond(cond)
                    

    #---------------------------------------------------------------------------          
    # Delete a projector
    #---------------------------------------------------------------------------      
    def __del__(self):
        qoptkit.prj_destroy_projector(c_long(self.obj))

    #---------------------------------------------------------------------------          
    # Set the post-selection condition from a vector
    #---------------------------------------------------------------------------      
    def set_cond(self, vec):   
        length=len(vec)
        send=(c_int*length)()
        for i in range(0,length):
            send[i]=vec[i]
            
        qoptkit.prj_add_term(c_long(self.obj),c_double(1.0), c_double(0.0), send) 
        
        
#------------------------------------------------------------------------------#
# Quantum optical circuit simulator                                            #
#------------------------------------------------------------------------------#     
class simulator(object):
    """

    Simulator that can be used to calculate the output state of a circuit given an input state.

    :mem(optional([int]): Reserved memory for the output expressed as a maximum number of kets. (Internal memory)
    
    """
    #---------------------------------------------------------------------------          
    # Create a simulator
    #---------------------------------------------------------------------------      
    def __init__(self, mem=1000):
        func=qoptkit.sim_new_simulator
        func.restype=c_long    
        self.obj = func(mem)

    #---------------------------------------------------------------------------      
    # Delete a simulator
    #---------------------------------------------------------------------------      
    def __del__(self):
        qoptkit.sim_destroy_simulator(c_long(self.obj))

    #---------------------------------------------------------------------------      
    # Run a simulator for an input state and a circuit
    #---------------------------------------------------------------------------      
    def run(self, istate,qoc):
        """

        Calculates an output state as a function of an input initial state for a given quantum circuit.

        :istate(state): Initial state.
        :qoc(qocircuit): Input quantum circuit.
        :return(state): Output state.
    
        """
        func=qoptkit.sim_run
        func.restype=c_long
        obj=func(c_long(self.obj),c_long(istate.obj),c_long(qoc.obj)) 
        newstate=state(1,dummy=True)
        newstate.obj=obj
        return newstate


        
#------------------------------------------------------------------------------#              
# QOPTKIT-QISKIT Compiler                                                        #
#------------------------------------------------------------------------------#      
class compiler(object):
    """

    Compiler between optical circuits defined in qoptkit and logical circuits defined in QISKIT.


    :mem(optional([int]): Number of kets in the internal memory.
             
    """
    #---------------------------------------------------------------------------      
    # Create a compiler
    #---------------------------------------------------------------------------      
    def __init__(self, mem=1000):
        self.sim=simulator(mem)        
       
    #---------------------------------------------------------------------------      
    # Creates an unitary matrix in qubit encoding representative of the optical 
    # circuit functionality. Auxilary method. Intended for internal use.
    #---------------------------------------------------------------------------   
    def circuit_mat(self, in_cond, out_cond, occ, psel, qoc, norm):
        # Initialize parameterss
        param1=to_int_ptr(in_cond)
        param2=to_int_ptr(out_cond)
        ncoef=2**param1[2]


        length=len(occ)
        send1=(c_int*length)()
        for i in range(0,length):
            send1[i]=occ[i]

        length=len(psel)
        send2=(c_int*length)()
        for i in range(0,length):
            send2[i]=psel[i]

        
        # Call C++ function
        func=qoptkit.qk_assemble
        func.restype=POINTER(c_double)
        array_ptr=func(c_long(self.sim.obj),param1[0],param2[0],param1[1],send1,send2,c_long(qoc.obj),norm)
        
        # Copy data and free memory
        ArrayType = 2*ncoef*ncoef*c_double
        addr = ctypes.addressof(array_ptr.contents)
        aux = np.frombuffer(ArrayType.from_address(addr), dtype=complex)       
        aux.shape=(ncoef,ncoef)
        array=aux.copy() 
        free_ptr(array_ptr)

        # Return results
        
        return array, param1[2]
    

    #---------------------------------------------------------------------------      
    # Returns a non-transpiled circuit. Auxilary method. Intended for internal use.
    #---------------------------------------------------------------------------      
    def circuit(self, in_cond, out_cond, occ, psel, qoc):
        mtx, nq = self.circuit_mat(in_cond,out_cond,occ, psel, qoc,1)
        circ = QuantumCircuit(nq,nq)
        gate = UnitaryGate(mtx)        
        circ.append(gate, list(range(nq-1,-1,-1)))
               
        return circ, nq
    
    

    #---------------------------------------------------------------------------      
    # Assemble circuit. Auxilary method. Intended for internal use.
    #---------------------------------------------------------------------------      
    def assemble_circ(self,in_cond, out_cond, occ, psel, basis, qoc):
        circuit, nq=self.circuit(in_cond,out_cond, occ, psel, qoc)
        transpiled=transpile(circuit, basis_gates=basis)

        return transpiled, nq
    
    #---------------------------------------------------------------------------      
    # Assemble full circuit: Circuit + Measurement   
    #---------------------------------------------------------------------------      
    def assemble(self,in_cond, out_cond, occ, psel, basis, qoc ,measure="qoptkit"):
        """
        Creates a transpiled circuit geiven an optical circuit  |br|

        :in_cond(list[][]):  nx2 matrix with the input qubit definitions. 
        :out_cond(list[][]): nx2 matrix with the output qubit definitions.
        :occ(list[int]): List of default occupation by input channels in case the channels are not defined as part of in_cond.
        :psel(list[int]): Post-selection condition of the circuit as a list of occupations by channel. Negative occupation means no post-selection for that channel.
        :basis(list): List of the allowed gates for transpilation.
        :qoc(qocircuit): Circuit to be assembled
        :measure(optional[string]): We measure following "qoptkit" or "qiskit" protocols. 
                                    In qoptkit the smaller qubit (qubit 0) appear at the left  when defining the outcomes while in QISKIT is the other way around.
        :return circ(QuantumCircuit): Returns the full assembled circuit.
        
        """  
        
        # Append the main circuit
        compiled, nq= self.assemble_circ(in_cond,out_cond, occ, psel, basis, qoc)
                
        # Append measurement
        # We can choose little endian like in qoptkit or big endian style like in qiskit
        if measure=='qoptkit':
            compiled.barrier()
            compiled.measure(range(0,nq), range(nq-1,-1,-1))
        if measure=='qiskit':
            compiled.barrier()
            compiled.measure(range(0,nq), range(0,nq))
            
        return compiled
    
    #---------------------------------------------------------------------------      
    # Returns the circuit matrix
    #---------------------------------------------------------------------------      
    def return_matrix(self, in_cond, out_cond, occ, psel,qoc):
        """

         Returns the matrix that relates the qubit encoded input and output states to export to QISKIT. |br|


        :in_cond(list[][]):  nx2 matrix with the input qubit definitions. 
        :out_cond(list[][]): nx2 matrix with the output qubit definitions.
        :occ(list[int]): List of default occupation by input channels in case the channels are not defined as part of in_cond.
        :psel(list[int]): Post-selection condition of the circuit as a list of occupations by channel. Negative occupation means no post-selection for that channel.        
        :qoc(qocircuit): Circuit to obtain the logic matrix
        :return mtx(list[][]): Returns a square circuit matrix containing the matrix in qubit encoding that represents the logic of
                               the optical circuit. 
        
        """
        mtx,nq, ptr=self.circuit_mat(in_cond,out_cond,occ, psel, qoc,1)
        return mtx
        
    #---------------------------------------------------------------------------      
    # Transpile a logical circuit into a quantum circuit by substitution
    #---------------------------------------------------------------------------      
    def transpile(self, qc):
        """

         Transpile a logical circuit into a quantum circuit using a substitution method. |br|

        :qc(QuantumCircuit):  Circuit defined in QISKIT to be transpiled.
        :return (qocircuit, list[int][int] , list[int], QuantumCircuit): Returns in this order: |br|
                                                                    * An optical circuit fullfilling the same functionality that the logical circuit. |br|
                                                                    * A matrix with the qubit-channel encoding. |br|
                                                                    * A list of the ancilla channels. |br|
                                                                    * A list with the post-selection condition. |br|
                                                                    * A transpiled QuantumCircuit with rotations and controled phase flips as gates. |br|
        
        """
        # Transpile the circuit to gates implementable in optical circuits
        transpiled=transpile(qc, basis_gates={'rx','ry','rz','cz'})
        nqbits=qc.num_qubits
        
        # Pre-process
        # Create translation condition
        cond=np.zeros((2,nqbits),dtype=int) 
        ch=0
        for i in range(0,nqbits):
            cond[i][0]=ch
            cond[i][1]=ch+1
            ch=ch+2
            
        # Compute number of ancillas
        ncz=0
        for gate in transpiled.data:
            if gate[0].name=='cz':
                ncz=ncz+1
                
        # Process
        # Create circuit
        nch=2*nqbits+2*ncz
        fancilla=2*nqbits
        optical=qocircuit(nch)
        
        #Process circuit
        nancilla=0
        for gate in transpiled.data:
            if(gate[0].name=='rx'):      
                qubit=gate[1]
                idx=qubit[0].index
                theta=0.5*180*gate[0].params[0]/np.pi
                optical.bs(cond[idx][0],cond[idx][1],theta,180.0/2.0)

            if(gate[0].name=='ry'):      
                qubit=gate[1]
                idx=qubit[0].index
                theta=0.5*180*gate[0].params[0]/np.pi
                optical.bs(cond[idx][0],cond[idx][1],-theta,0.0)
            
            if(gate[0].name=='rz'):      
                qubit=gate[1]
                idx=qubit[0].index
                phi=180*gate[0].params[0]/np.pi
                optical.ph(cond[idx][0],  phi/2.0)
                optical.ph(cond[idx][1], -phi/2.0)

            if(gate[0].name=='cz'):  
                # Csign configuration
                theta1 = +54.74
                theta2 = -54.74
                theta3 = +54.74
                theta4 = +17.63
                phi1   = +0.0
                phi2   = +0.0
                phi3   = +0.0
                phi4   = +0.0
                
                # Qubits
                qubit=gate[1]
                q0=qubit[0].index
                q1=qubit[1].index
                
                # Translation to channels
                ch0=cond[q0][0]
                ch2=cond[q1][0]
                ch4=fancilla+nancilla;
                ch5=fancilla+nancilla+1                
                nancilla=nancilla+2                
                
                # CSIGN 2/27 circuit
                optical.barrier()
                optical.ph(ch0, 180.0)
                optical.ph(ch2, 180.0)
                
                optical.bs(ch0,ch4,theta1,phi1);
                optical.bs(ch2,ch5,theta3,phi3);
                
                optical.bs(ch0,ch2,theta2,phi2);
                optical.bs(ch4,ch5,theta4,phi4);
                optical.barrier()   
                

        
        # Set up detectors
        ancillas=[]
        pselect=[-1]*nch
        # Ancilla detectors
        for i in range(0,nancilla):
            ancillas.append(1)
            pselect[fancilla+i]=1
            
        # Return circuit (and other information)
        return optical, cond , ancillas, pselect, transpiled        
