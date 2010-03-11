#!/usr/bin/python
####
####
##
##  Original logic Andreas (thanks), rewritten in OOP, tested etc. by Leo
## 
## Copyright (C) 2002 Andreas Heger All rights reserved
## Copyright (C) 2007 Leo Goodstadt All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
## Author: Leo Goodstadt <leo.goodstadt@dpag.ac.uk>
##
##
##
####
####

import os,sys,time

# add self to search path for testing
if __name__ == '__main__':
    sys.path.append("../python_modules")

import general_util

from collections import namedtuple
t_interval = namedtuple("t_interval", "beg end")



#----------------------------------------------------------------
class t_intervals (object):
    def __init__ (self, interval_list=None):
        if interval_list == None:
            self.data = list()
        else:
            self.data = interval_list

    def __repr__ (self):
        return repr(self.data)
            
    def __cmp__ (self, other):
        return cmp(self.data, other.data)


    #----------------------------------------
    #   Exclusive 
    # 
    #       e.g. 1-10  11-12
    #       
    #       intervals are 2 units apart
    #       
    #       i.e. this_from - last_to + 1 = distance
    # 
    #            gap = (last_to, this_from)
    #            
    #   distance = this_from - last_to + 1
    # 
    #   overlap  = this_from  >  last_to - 1
    # 
    #   gap = (last_to, this_from)
    # 
    # 
    def combine_close_neighbours( self, min_distance ):
        """
            combine_close_neighbours
            
                Merge neighbours that are less than a certain
                distance apart.
        
                side effect: sort intervals
        
        """
        
        if not self.data or not len(self.data): return t_intervals([])
    
        new_data = []
        
        self.data.sort()
        
        first_from, last_to = self.data[0][0:2]
        extras = [self.data[0][2:]]
        has_extras = len(self.data[0]) > 2
            


        
        for d in self.data[1:]:
            this_from, this_to = d[0:2]

            if this_from - last_to + 1 > min_distance:
                if has_extras:
                    new_data.append(type(self.data[0])((first_from, last_to, extras) ))
                else:
                    new_data.append(type(self.data[0])((first_from, last_to) ))
                first_from = this_from
                extras = []
            else:
                extras.append(d[2:])
            last_to = max(this_to, last_to)
        
        if has_extras:
            new_data.append(type(self.data[0])((first_from, last_to, extras) ))
        else:
            new_data.append(type(self.data[0])((first_from, last_to) ))
    
        self.data = new_data

        return self


    #----------------------------------------
    def combine_overlapping (self):
        """
            combine_overlapping
    
                Overlapping intervals are concatenated into larger intervals.
    
                side effect: sort intervals
        """
        self.combine_close_neighbours(0 )

        return self
    

    #----------------------------------------
    def union (self, other):
        """
            combine_overlapping
    
                Overlapping intervals are concatenated into larger intervals.
    
                side effect: sort intervals
                """
        return t_intervals(self.data + other.data).combine_overlapping()
    



    #----------------------------------------
    def complemented(self, first = None, last = None):
        """
            complemented
        
                returns sets of intervals between gaps of
                current intervals
        
                if first and last are supplied,
                will union opening and ending gaps as well
        
                side effect: sort intervals
            
                assumes combined
        """
    
        if not self.data:
            if first and last:
                return t_intervals([t_interval(first, last)])
            else:
                return t_intervals([])
        
        new_data = []
        
        self.data.sort()
        last_from, last_to = self.data[0][0:2]
    
        if first != None and last_from  >  first - 1:
            new_data.append( t_interval(first, last_from) )
    
        for d in self.data:
            this_from, this_to = d[0:2]
            if this_from > last_to - 1:
                new_data.append( t_interval(last_to  , 
                                  this_from ) )            
    
            last_from = this_from
            last_to = max(last_to, this_to)
    
        if last != None:
            if last  >  last_to - 1:
                new_data.append( t_interval(last_to, last))
    
        return t_intervals(new_data)


    #----------------------------------------
    def subtract( self, to_remove ):
        """
        remove intervals contained in to_remove
    
        assumes combined
        
        """
        if not to_remove.data or not self.data:
            return t_intervals([])

        
        new_data = []
        
        self.data.sort()
        to_remove.data.sort()
    
        current_to_remove = 0
        
        for d in self.data:
            this_from, this_to = d[0:2]
    
            for dd in to_remove.data:
                remove_from, remove_to = dd[0:2]
                #debug
                #print (this_from, this_to), (remove_from, remove_to)

                #this_from  >  last_to - 1
                #last_to - 1 < this_from

                # remove_to does not overlap with this interval
                if remove_to - 1 < this_from: continue
                if remove_from > this_to - 1: continue
    

                # contained entirely within: kill this interval
                if remove_from <= this_from and remove_to >= this_to:
                    this_from = remove_to # mark for death
                    break

                # loped off end of this interval
                if this_from < remove_from:
                    new_data.append( t_interval(this_from, remove_from) )
                    # print "unioning", this_from, remove_from
                
                this_from = max(this_from, remove_to)
    
                if this_to - 1 <= this_from: break
                
            if this_to - 1 > this_from:
                # print "unioning", this_from, this_to
                new_data.append( t_interval(this_from, this_to) )            
    
    
        return t_intervals(new_data)




    

    #----------------------------------------
    def with_small_intervals_filtered(self, min_length):
        """
        with_small_intervals_filtered
    
            delete intervals shorter than a minimum length
        """
        
        if not self.data: return t_intervals([])
    
        return t_intervals([d for d in self.data if d[1] - d[0] >= min_length])


    #----------------------------------------------------------------
    #----------------------------------------------------------------
    def calculate_overlap_with( self, other ):
        """
        calculate_overlap_with

            assumes combined

            returns overlap
            
        """
    
        if not self.data or not other.data:
            return 0
    
        self.data.sort()
        other.data.sort()
    
        overlap = 0
        x = 0
        y = 0
    
        while x < len(self.data) and y < len(other.data):
    
            xfrom, xto = self.data[x] [0:2]
            yfrom, yto = other.data[y][0:2]
    
            if xto - 1 < yfrom:
                x += 1
            elif yto - 1 < xfrom:
                y += 1
            else:
                overlap += min( xto, yto ) - max(xfrom, yfrom )
                
                if xto < yto:
                    x += 1
                elif yto < xto:
                    y += 1
                else:
                    x += 1
                    y += 1
    
        return overlap

    #----------------------------------------------------------------
    #----------------------------------------------------------------
    def intersect( self, other ):
        """
        intersect

            assumes combined

            returns overlap
            
        """
    
        if not self.data or not other.data:
            return t_intervals([])
    
        self.data.sort()
        other.data.sort()
    
        overlap = 0
        x = 0
        y = 0

        new_data = []
        while x < len(self.data) and y < len(other.data):
    
            xfrom, xto = self .data[x][0:2]
            yfrom, yto = other.data[y][0:2]
    
            if xto - 1 < yfrom:
                x += 1
            elif yto - 1 < xfrom:
                y += 1
            else:
                new_data.append((max(xfrom, yfrom), min( xto, yto )))
                
                if xto < yto:
                    x += 1
                elif yto < xto:
                    y += 1
                else:
                    x += 1
                    y += 1
    
        return t_intervals(new_data)

    #----------------------------------------------------------------
    #----------------------------------------------------------------
    def total_len( self):
        """
        total_len

            assumes combined

            cummulative lengths
            
        """
    
        if not self.data:
            return 0
    
        return sum(d[1] - d[0] for d in self.data)
    
    #----------------------------------------------------------------
    #----------------------------------------------------------------
    def find_containing_interval( self, test_pt):
        """
        find_interval

            assumes combined and sorted

            finds first interval which contains point
            
        """
        #print test_pt
        if not self.data:
            return (None, None)
    
        for d in self.data:
            fro, to = d[0:2]
            if test_pt < fro:
                #print "break"
                break
            if test_pt > to - 1:
                #print "next"
                continue
            #print "yes"
            return (fro, to)
            
        return (None, None)
    
 

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


#        unit test
    

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
def test_combine ():
    # test for non-combine of exclusive
    intervals1 = t_intervals([(1,10), (10, 40), (30, 70) ])
    intervals1.combine_overlapping()
    assert(intervals1 == t_intervals([(1,10), (10, 70)]))

    # test for combine of exclusive
    intervals1 = t_intervals([(1,10), (9, 40), (30, 70) ])
    intervals1.combine_overlapping()
    assert(intervals1 == t_intervals([(1, 70)]))



def test_combine_close_neighbours ():
    # test for non-combine of exclusive
    intervals1 = t_intervals([(1,10), (12, 40), (30, 70) ])
    intervals1.combine_close_neighbours(2)
    assert(intervals1 == t_intervals([(1,10), (12, 70)]))

    # test for combine of exclusive
    intervals1 = t_intervals([(1,10), (11, 40), (30, 70) ])
    intervals1.combine_close_neighbours(2)
    assert(intervals1 == t_intervals([(1, 70)]))


def test_complement ():
    # test exclusive complement
    intervals1 = t_intervals([(4,10), (12, 40), (50, 70) ])
    assert(intervals1.complemented() == t_intervals([(10, 12), (40, 50)]) )

    # test exclusive complement with ranges
    intervals1 = t_intervals([(4,10), (12, 40), (50, 70) ])
    assert(intervals1.complemented(1, 200) == t_intervals([(1, 4), (10, 12), (40, 50), (70, 200)]) )


def test_subtract ():
    # test exclusive
    intervals1 = t_intervals([(1,20), (50, 100) ])
    intervals2 = t_intervals([(15,55), (70, 85) ])
    assert intervals1.subtract(intervals2) == t_intervals([(1, 15), (55, 70), (85, 100)])
    

def test_with_small_intervals_filtered ():
    # test exclusive
    intervals1 = t_intervals([(1,10), (20,30),  (40, 51) ])
    assert(intervals1.with_small_intervals_filtered(11) == t_intervals( [(40, 51)]))

           

def test_overlap ():
    # test exclusive complement
    intervals1 = t_intervals([(1,20), (50, 100) ])
    intervals2 = t_intervals([(15,55), (70, 85) ])
    assert(intervals1.calculate_overlap_with(intervals2) == 25)
    assert intervals1.intersect(intervals2).total_len() ==intervals1.calculate_overlap_with(intervals2)

def test_total_len():
    intervals1 = t_intervals([(1,20), (50, 100) ])
    intervals2 = t_intervals([(15,55), (70, 85) ])
    assert (intervals1.total_len() + 
           intervals2.total_len() -
           intervals1.calculate_overlap_with(intervals2) == 
           intervals1.union(intervals2).total_len())


def test_find_containing_interval():
    intervals1 = t_intervals([(1,20), (50, 100) , (140, 160)])
    intervals1.combine_overlapping()
    assert (intervals1.find_containing_interval(0) == (None,None)) 
    assert (intervals1.find_containing_interval(1) == (1,20)) 
    assert (intervals1.find_containing_interval(19) == (1,20)) 
    assert (intervals1.find_containing_interval(20) == (None, None)) 

    assert (intervals1.find_containing_interval(49) == (None, None)) 
    assert (intervals1.find_containing_interval(50) == (50, 100)) 
    assert (intervals1.find_containing_interval(99) == (50, 100)) 
    assert (intervals1.find_containing_interval(100) == (None, None)) 

    assert (intervals1.find_containing_interval(139) == (None, None)) 
    assert (intervals1.find_containing_interval(140) == (140, 160)) 
    assert (intervals1.find_containing_interval(159) == (140, 160)) 
    assert (intervals1.find_containing_interval(160) == (None, None)) 


def test_intersect():
    intervals1 = t_intervals([(1,20), (50, 100) ])
    intervals2 = t_intervals([(15,55), (70, 85) ])
    assert (intervals1.intersect(intervals2) == t_intervals([(15, 20), (50, 55), (70, 85)]))


def unit_test ():
    test_combine()
    test_combine_close_neighbours ()
    test_complement ()
    test_subtract()
    test_with_small_intervals_filtered ()
    test_overlap()
    test_total_len()
    test_find_containing_interval()
    test_intersect()
    print "All tests pass"

if __name__ == '__main__':
    unit_test()

