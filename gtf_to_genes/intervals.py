#!/usr/bin/python
####
####
##  Latest code
##
##  Some of original logic Andreas (thanks), rewritten in OOP, tested etc. by Leo
##
## Copyright (C) 2002 Andreas Heger All rights reserved
## Copyright (C) 2007-2014 Leo Goodstadt All rights reserved
##
## Author: Andreas Heger <heger@dpag.ox.ac.uk>
## Author: Leo Goodstadt <python@llew.org.uk>
##
## Added support for additional fields
#       N.B. Use tuples internally
##
####
####

import os,sys,time

from collections import namedtuple

#
#   only used in complemented
#
t_interval = namedtuple("t_interval", "beg end")


def make_interval(begin, end, curr_interval, istuple, data_type):
    if istuple:
        return (begin, end) + curr_interval[2:]
    else:
        return data_type(*(begin, end) + curr_interval[2:])


#----------------------------------------------------------------
class t_intervals (object):
    def append (self, item):
        if not isinstance(item, tuple):
            raise Exception("%s is not a tuple" % str(item))
        self.data.append(item)
    def __init__ (self, interval_list=None):
        if interval_list is None:
            self.data = list()
        else:
            self.data = interval_list

    def __repr__ (self):
        return repr(self.data)

    def __cmp__ (self, other):
        return cmp(self.data, other.data)


    #----------------------------------------
    def remove_contained_overlaps(self):
        """
            remove intervals that are fully contained in another.

            [(10, 100), (20, 50), (70, 120), (130, 200), (10, 50), (140, 210), (150, 200)]

            results:

            [(10, 100), (70, 120), (130, 200), (140, 210)]

        """
        if not self.data or not len(self.data):
            return self

        new_data = []

        self.data.sort()

        prev_index = 0

        for curr_index in range(1, len(self.data)):

            # prev is larger, drop curr_index
            if self.data[prev_index][1] >= self.data[curr_index][1]:
                continue

            # current is longer and completely overlaps previous, drop prev_index
            if self.data[prev_index][0] == self.data[curr_index][0]:
                prev_index = curr_index
                continue

            # no complete overlap
            new_data.append( self.data[prev_index])

            prev_index = curr_index

        new_data.append( self.data[prev_index])

        self.data = new_data

        return self

    #----------------------------------------
    def shorten_at_overlap(self):
        """
            1) distribute overlaps between neighbouring intervals
            2) remove intervals that are fully contained in another.

            N.B. This is a greedy algorithm and may give consistent but unexpected
            results when three successive intervals overlap...

        """
        if not self.data or not len(self.data):
            return self

        new_data = []

        self.data.sort()
        #print >>sys.stderr, self.data

        prev = self.data[0]

        minimum = self.data[0][0]

        # to recreate datatype
        data_type = type(self.data[0])
        istuple = isinstance(self.data[0], tuple)

        for curr in self.data[1:]:

            # curr should be limited by previous guillotine
            if curr[0] < minimum:
                # discard curr altogether
                if curr[1] < minimum:
                    #print >>sys.stderr, "Drop curr %s < minimum (%d)" % (curr, minimum)
                    continue
                #print >>sys.stderr, "Adjust curr %s to minimum (%d)" % (curr, minimum)
                curr = make_interval(minimum, curr[1], curr[2:], istuple, data_type)


            # prev is larger, drop curr
            if prev[1] >= curr[1]:
                #print >>sys.stderr, "prev %s is larger, drop curr %s" % (prev, curr)
                continue

            # current is longer and completely overlaps previous, drop prev
            if prev[0] == curr[0]:
                #print >>sys.stderr, "curr %s is longer and completely overlaps prev, drop prev %s" % (curr, prev)
                prev = curr
                continue

            # if prev and curr do not overlap, add prev
            if prev[1] <= curr[0]:
                new_data.append(prev)
                minimum = curr[0]
                #print >>sys.stderr, "Curr %s does not overlap prev %s, add prev, minimum = %s" % (curr, prev, minimum)
                prev = curr
                continue

            # overlap: shorten both
            midpoint = int((prev[1] - curr[0]) / 2) + curr[0]
            #print >>sys.stderr, "overlap: shorten curr %s, prev %s to midpoint %s " % (curr, prev, midpoint)
            prev = make_interval(prev[0], midpoint, prev, istuple, data_type)
            curr = make_interval(midpoint, curr[1], curr, istuple, data_type)
            new_data.append(prev)
            prev = curr
            minimum = midpoint

        new_data.append(prev)

        self.data = new_data

        return self


    #----------------------------------------
    def extend_into_gaps(self, minimum, maximum):
        """
            Extend neighbouring intervals into gaps

        """
        if not self.data or not len(self.data):
            return self

        new_data = []

        self.data.sort()
        #print >>sys.stderr, self.data

        # to recreate datatype
        data_type = type(self.data[0])
        istuple = isinstance(self.data[0], tuple)

        prev_max = self.data[0]

        if prev_max[0] > minimum:
            prev_max = make_interval(minimum, prev_max[1], prev_max, istuple, data_type)

        for curr in self.data[1:]:

            # curr touches prev
            if curr[0] < prev_max[1]:
                # curr only overlaps prev
                if curr[1] > prev_max[1]:
                    new_data.append(prev_max)
                    prev_max = curr
                # prev encompasses curr
                else:
                    new_data.append(curr)
                    # keep prev_max
                continue

            # doesn't touch:
            midpoint = int((curr[0] + prev_max[1]) / 2)
            # extend both prev and append
            new_data.append(make_interval(prev_max[0], midpoint, prev_max, istuple, data_type))
            # extend curr which becomes the new prev_max
            prev_max = make_interval(midpoint, curr[1], curr, istuple, data_type)

        if prev_max[1] < maximum:
            new_data.append(make_interval(prev_max[0], maximum, prev_max, istuple, data_type))
        else:
            new_data.append(prev_max)
        new_data.sort()

        self.data = new_data

        return self


    #----------------------------------------
    #   open-closed notation
    #
    #       e.g. [1-10)  [11-12)
    #
    #       intervals are 1 unit apart
    #
    #       i.e. this_beg - last_end = distance
    #
    #            gap = (last_end, this_beg)
    #
    #   distance = this_beg - last_end
    #
    #   overlap  = this_beg  >=  last_end
    #
    #   gap = [last_end, this_beg)
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

        data_type = type(self.data[0])
        istuple = isinstance(self.data[0], tuple)

        curr_beg, curr_end = self.data[0][0:2]
        extras = self.data[0]

        #print >>sys.stderr, "1", (curr_beg, curr_end)
        for d in self.data[1:]:
            next_beg, next_end = d[0:2]
            #print >>sys.stderr, "-", (next_beg, next_end)

            # Non - overlap
            if next_beg - curr_end >= min_distance:
                #print >>sys.stderr, "+", (curr_beg, curr_end), (next_beg, next_end)
                new_data.append(make_interval(curr_beg, curr_end, extras, istuple, data_type))
                curr_beg = next_beg
                extras = d

            curr_end = max(next_end, curr_end)

        new_data.append(make_interval(curr_beg, curr_end, extras, istuple, data_type))

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
    def remove_empty (self):
        """
            remove_empty

                Remove empty intervals

                side effect: sort intervals
        """
        self.data = [d for d in self.data if d[1] > d[0]]


        return self

    #----------------------------------------
    def unioned_with (self, other):
        """
            combine_overlapping

                Overlapping intervals are concatenated into larger intervals.

                side effect: sort intervals
                """
        return t_intervals(self.data + other.data).combine_overlapping()

    def union( self, other ):
        new_data = self.unioned_with(other)
        self.data = new_data.data
        return self

    def __or__(self, other) :
        """
        self.union(other)
        """
        return self.union(other)

    def __ior__(self, other) :
        """
        self.unioned_with(other)
        """
        return self.unioned_with (other)


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
        last_beg, last_end = self.data[0][0:2]

        if first != None and last_beg  >  first - 1:
            new_data.append( t_interval(first, last_beg) )

        for d in self.data:
            this_beg, this_end = d[0:2]
            if this_beg > last_end - 1:
                new_data.append( t_interval(last_end, this_beg ) )

            last_beg = this_beg
            last_end = max(last_end, this_end)

        if last != None:
            if last  >  last_end - 1:
                new_data.append( t_interval(last_end, last))

        return t_intervals(new_data)

    #----------------------------------------
    def expanded_by (self, size):
        """
            combine_overlapping

                Overlapping intervals are concatenated into larger intervals.

                side effect: sort intervals
        """
        if not self.data or not len(self.data): return t_intervals([])

        new_data = []
        data_type = type(self.data[0])
        istuple = isinstance(self.data[0], tuple)

        for d in self.data:
            new_beg, new_end = d[0:2]
            new_beg = max(0,new_beg - size)
            new_end += size

            new_data.append(make_interval(new_beg, new_end, d, istuple, data_type))

        return t_intervals(new_data).combine_overlapping()

    def expand_by( self, size ):
        new_data = self.expand_by(size)
        self.data = new_data.data
        return self

    def __iadd__(self, size):
        """
        self.expanded_by(size)
        """
        return self.expand_by(size)

    def __add__(self, size):
        """
        self.expand_by(size)
        """
        return self.expanded_by(size)

    #----------------------------------------
    def subtracted_by( self, to_remove ):
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

        data_type = type(self.data[0])
        istuple = isinstance(self.data[0], tuple)
        for d in self.data:
            this_beg, this_end = d[0:2]

            for dd in to_remove.data:
                remove_beg, remove_end = dd[0:2]
                #debug
                #print (this_beg, this_end), (remove_beg, remove_end)

                #this_beg  >  last_end - 1
                #last_end - 1 < this_beg

                # remove_end does not overlap with this interval
                if remove_end - 1 < this_beg: continue
                if remove_beg > this_end - 1: continue


                # contained entirely within: kill this interval
                if remove_beg <= this_beg and remove_end >= this_end:
                    this_beg = remove_end # mark for death
                    break

                # loped off end of this interval
                if this_beg < remove_beg:
                    new_data.append(make_interval(this_beg, remove_beg, d, istuple, data_type))


                    # print "unioning", this_beg, remove_beg

                this_beg = max(this_beg, remove_end)

                if this_end - 1 <= this_beg: break

            if this_end - 1 > this_beg:
                # print "unioning", this_beg, this_end
                new_data.append(make_interval(this_beg, this_end, d, istuple, data_type))



        return t_intervals(new_data)


    def subtract( self, other ):
        new_data = self.subtracted_by(other)
        self.data = new_data.data
        return self

    def __isub__(self, other):
        """
        self.subtracted_by(other)
        """
        return self.subtract(other)

    def __sub__(self, other):
        """
        self.subtract(other)
        """
        return self.subtracted_by(other)




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

        #print >> sys.stderr, self.data, other.data
        while x < len(self.data) and y < len(other.data):

            xbeg, xend = self.data[x] [0:2]
            ybeg, yend = other.data[y][0:2]

            if xend - 1 < ybeg:
                x += 1
            elif yend - 1 < xbeg:
                y += 1
            else:
                overlap += min( xend, yend ) - max(xbeg, ybeg )

                if xend < yend:
                    x += 1
                elif yend < xend:
                    y += 1
                else:
                    x += 1
                    y += 1

        return overlap

    #----------------------------------------------------------------
    #----------------------------------------------------------------
    def intersected_with( self, other ):
        """
        intersect

            assumes combined

            returns overlap

        """

        if not self.data or not other.data:
            return t_intervals([])

        data_type = type(self.data[0])
        istuple = isinstance(self.data[0], tuple)


        self.data.sort()
        other.data.sort()

        overlap = 0
        x = 0
        y = 0

        new_data = []
        while x < len(self.data) and y < len(other.data):

            xbeg, xend = self .data[x][0:2]
            ybeg, yend = other.data[y][0:2]

            if xend - 1 < ybeg:
                x += 1
            elif yend - 1 < xbeg:
                y += 1
            else:
                # remember to add self. extra fields
                new_data.append(make_interval(max(xbeg, ybeg), min( xend, yend ), self.data[x], istuple, data_type))

                if xend < yend:
                    x += 1
                elif yend < xend:
                    y += 1
                else:
                    x += 1
                    y += 1

        return t_intervals(new_data)

    def intersect( self, other ):
        new_data = self.intersected_with(other)
        self.data = new_data.data
        return self

    def __iand__(self, other):
        """
        self.intersected_with(other)
        """
        return self.intersect(other)

    def __and__(self, other):
        """
        self.intersect(other)
        """
        return self.intersected_with(other)



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
            beg, end = d[0:2]
            if test_pt < beg:
                #print "break"
                break
            if test_pt > end - 1:
                #print "next"
                continue
            #print "yes"
            return d

        return (None, None)

