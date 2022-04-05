#!/usr/bin/env python
# coding: utf-8

# # Two ways of implementing a 3-state cellular automata with the neighborhood defined as the state cell directly above and the cell above and to the left

# ## (adapted from cellular automata code created by Adam Rupe)

# ## Imports

# In[1]:


import random
from matplotlib import pyplot as plt


# ## Spaghetti code implementation of 3-state cellular automata (no functions or classes)

# In[3]:


rule_number = 8711
length = 100
time = 100

# make the initial condition
initial_condition = []
for i in range(length):
    initial_condition.append(random.randint(0,2))

# create list of neighborhood tuples in lex. order
neighborhoods = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)]

# convert the rule number to ternary and pad with 0s as needed
tern_nums = []
while rule_number > 0:
    rule_number, r = divmod(rule_number,3)
    tern_nums.append(str(r))
in_ternary = ''.join(tern_nums)
    
ternary_length = len(in_ternary)
if ternary_length != 9:
    padding = 9 - ternary_length
    in_ternary = in_ternary + '0'*padding

# create the lookup table dictionary
lookup_table = {}
for i in range(9):
    key = neighborhoods[i]
    val = in_ternary[i]
    lookup_table.update({key:val})
    
# initialize spacetime field and current configuration
spacetime_field = [initial_condition]
current_configuration = initial_condition.copy()

# apply the lookup table to evolve the CA for the given number of time steps
for t in range(time):
    new_configuration = []
    for i in range(len(current_configuration)):
        
        neighborhood = (current_configuration[(i-1)], 
                        current_configuration[i])
        
        new_configuration.append(int(lookup_table[neighborhood]))
        
    current_configuration = new_configuration
    spacetime_field.append(new_configuration)
    
# plot the spacetime field diagram
plt.figure(figsize=(12,12))
plt.imshow(spacetime_field, cmap=plt.cm.Greys, interpolation='nearest')
plt.show()


# ## Implementation of 3-state cellular automata using functions and classes

# In[4]:


def random_ternary_string(length):
    '''
    Returns a random ternary string of the given length. 
    
    Parameters
    ----------
    length: int
        Posivite integer that specifies the desired length of the ternary string.
        
    Returns
    -------
    out: list
        The random ternary string given as a list, with int elements.
    '''
    if not isinstance(length, int) or length < 0:
        raise ValueError("input length must be a positive ingeter")
    return [random.randint(0,2) for _ in range(length)]


# In[5]:


def three_state_lookup_table(rule_number):
    '''
    Returns a dictionary which maps 3-state cellular automata neighborhoods to output values. 
    Uses Wolfram rule number convention.
    
    Parameters
    ----------
    rule_number: int
        Integer value between 0 and 19682, inclusive. Specifies the CA lookup table
        according to the Wolfram numbering scheme.
        
    Returns
    -------
    lookup_table: dict
        Lookup table dictionary that maps neighborhood tuples to their output according to the 
        ECA local evolution rule (i.e. the lookup table), as specified by the rule number. 
    '''
    if not isinstance(rule_number, int) or rule_number < 0 or rule_number > 19682:
        raise ValueError("rule_number must be an int between 0 and 19682, inclusive")
    neighborhoods = [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)]
    
    tern_nums = []
    while rule_number > 0:
        rule_number, r = divmod(rule_number,3)
        tern_nums.append(str(r))
    tern_num = ''.join(reversed(tern_nums))
    
    in_ternary = '{:{fill}{align}{width}}'.format(tern_num, 
                                                  fill='0', 
                                                  align='>', 
                                                  width='9')
    
    return dict(zip(neighborhoods, map(int,reversed(in_ternary)))) # use map so that outputs are ints, not strings


# In[6]:


def three_state_spacetime_field(rule_number, initial_condition, time_steps):
    '''
    Returns a spacetime field array using the given rule number on the 
    given initial condition for the given number of time steps.
    
    Parameters
    ----------
    rule_number: int
        Integer value between 0 and 19682, inclusive. Specifies the three-state CA lookup table
        according to the Wolfram numbering scheme.
    initial_condition: list
        Binary string used as the initial condition for the three-state CA. Elements of the list
        should be ints. 
    time_steps: int
        Positive integer specifying the number of time steps for evolving the ECA. 
    '''
    if time_steps < 0:
        raise ValueError("time_steps must be a non-negative integer")
    # try converting time_steps to int and raise a custom error if this can't be done
    try:
        time_steps = int(time_steps)
    except ValueError:
        raise ValueError("time_steps must be a non-negative integer")
        
    # we will see a cleaner and more efficient way to do the following when we introduce numpy
    for i in initial_condition:
        if i not in [0,1,2]:
            raise ValueError("initial condition must be a list of 0s, 1s and 2s")
        
    lookup = three_state_lookup_table(rule_number)
    length = len(initial_condition)
    
    # initialize spacetime field and current configuration
    spacetime_field = [initial_condition]
    current_configuration = initial_condition.copy()

    # apply the lookup table to evolve the CA for the given number of time steps
    for t in range(time_steps):
        new_configuration = []
        for i in range(length):

            neighborhood = (current_configuration[(i-1)], 
                            current_configuration[i])

            new_configuration.append(lookup[neighborhood])

        current_configuration = new_configuration
        spacetime_field.append(new_configuration)
    
    return spacetime_field


# In[7]:


def spacetime_diagram(spacetime_field, size=12, colors=plt.cm.Greys):
    '''
    Produces a simple spacetime diagram image using matplotlib imshow with 'nearest' interpolation.
    
   Parameters
    ---------
    spacetime_field: array-like (2D)
        1+1 dimensional spacetime field, given as a 2D array or list of lists. Time should be dimension 0;
        so that spacetime_field[t] is the spatial configuration at time t. 
        
    size: int, optional (default=12)
        Sets the size of the figure: figsize=(size,size)
    colors: matplotlib colormap, optional (default=plt.cm.Greys)
        See https://matplotlib.org/tutorials/colors/colormaps.html for colormap choices.
        A colormap 'cmap' is called as: colors=plt.cm.cmap
    '''
    plt.figure(figsize=(size,size))
    plt.imshow(spacetime_field, cmap=colors, interpolation='nearest')
    plt.show()


# In[8]:


class three_state_CA(object):
    '''
    Three-state cellular automata simulator.
    '''
    def __init__(self, rule_number, initial_condition):
        '''
        Initializes the simulator for the given rule number and initial condition.
        
        Parameters
        ----------
        rule_number: int
            Integer value between 0 and 19682, inclusive. Specifies the CA lookup table
            according to the Wolfram numbering scheme.
        initial_condition: list
            Ternary string used as the initial condition for the three-state CA. Elements of the list
            should be ints. 
        
        Attributes
        ----------
        lookup_table: dict
            Lookup table for the three-state CA given as a dictionary, with neighborhood tuple keys. 
        initial: array_like
            Copy of the initial conditions used to instantiate the simulator
        spacetime: array_like
            2D array (list of lists) of the spacetime field created by the simulator.
        current_configuration: array_like
            List of the spatial configuration of the ECA at the current time
        '''
        # we will see a cleaner and more efficient way to do the following when we introduce numpy
        for i in initial_condition:
            if i not in [0,1,2]:
                raise ValueError("initial condition must be a list of 0s, 1s, and 2s")
                
        self.lookup_table = three_state_lookup_table(rule_number)
        self.initial = initial_condition
        self.spacetime = [initial_condition]
        self.current_configuration = initial_condition.copy()
        self._length = len(initial_condition)

    def evolve(self, time_steps):
        '''
        Evolves the current configuration of the three-state CA for the given number of time steps.
        
        Parameters
        ----------
        time_steps: int
            Positive integer specifying the number of time steps for evolving the three-state CA.  
        '''
        if time_steps < 0:
            raise ValueError("time_steps must be a non-negative integer")
        # try converting time_steps to int and raise a custom error if this can't be done
        try:
            time_steps = int(time_steps)
        except ValueError:
            raise ValueError("time_steps must be a non-negative integer")

        for _ in range(time_steps): # use underscore if the index will not be used
            new_configuration = []
            for i in range(self._length):

                neighborhood = (self.current_configuration[(i-1)], 
                                self.current_configuration[i])

                new_configuration.append(self.lookup_table[neighborhood])

            self.current_configuration = new_configuration
            self.spacetime.append(new_configuration)


# In[9]:


CA = three_state_CA(8711,random_ternary_string(100))
CA.evolve(100)
spacetime_diagram(CA.spacetime,10)


# In[ ]:




