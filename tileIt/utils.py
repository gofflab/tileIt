#!/usr/bin/env python
import sys,types,string

###################################
#
#Boolean Functions
#
####################################
def ifab(test, a, b):
    """x = ifab(test, a, b)
       WARNING:  Both 'a' and 'b' are evaluated
       C equivalent: x = test?a:b;
       Scheme equiv: (set x (if test a b))
       Python equiv: test and a or b
       None of the equivalents evaluates both arguments
    """
    if test: return a
    else: return b

###################################
#
#String Functions
#
####################################
def sfill(s, length, fill_char = '.'):
    #  Appends fill_char to the string s until it reaches length length
    #  ex:  sfill('hello',18,'.') -> hello...............
    #                                <---  18 chars  --->
    # useful for printing dictionaries in a cute way
    #    one......: 1
    #    five.....: 5
    #    seventeen: 17


    #list = map(None, s)
    #list.extend(map(None, fill_char*(length - len(list))))
    #return string.join(list, '')

    return s + fill_char*(length-len(s))


########
#
#Pretty Printing
#
########
def pretty_print(f, d, level=-1, maxw=0, maxh=0, gap="", first_gap='', last_gap=''):
    # depending on the type of expression, it recurses through its elements
    # and prints with appropriate indentation

    # f   is the output file stream
    # d   is the data structure
    #
    # level is the number of allowed recursive calls, the depth at which
    #       the data structure is explored
    #       default: -1 means never stop recursing early
    # maxw  is the maximum width that will be printed from the last element
    #       of the recursion (when no further recursion is possible, or
    #       the maximal depth has been reached)
    #       default: 0 means every line will be printed in its entirety, regardless
    #                of how long it may be
    # maxh  (max height) is the maximum number of elements that will be
    #       printed from a list or a dictionary, at any level or recursion
    #       default: 0 means every list or dictionary will have all its elements
    #                printed, even if it contains thousands of elements
    #
    # gap is the gap to include before every element of a list/dic/tuple
    # first_gap is the opening gap before the opening bracket, parens or curly braces
    # first_gap is the closing gap before the closing bracket, parens or curly braces
    
    if level == 0:
        if type(d) != types.StringType: d = `d`

        if maxw and len(d) > maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'...'+d[-final:]+' (%s chars)\n' % len(d))
        else: f.write(first_gap+d+'\n')
    elif type(d) == types.ListType:
        if not d:
            f.write(first_gap+"[]\n")
            return
        # recurse on lists
        f.write(first_gap+"[\n")
        h = 0
        for el in d:
            pretty_print(f, el, level-1, maxw, maxh, gap+'   ', gap+' ->', gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' -> ... (%s in list)\n'%len(d))
                    break
        f.write(last_gap+"]\n")
    elif type(d) == types.TupleType:
        if not d:
            f.write(first_gap+"()\n")
            return
        # recurse on tuples
        f.write(first_gap+"(\n")
        h = 0
        for el in d:
            pretty_print(f, el,
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'   ',
                         first_gap = gap+' =>',
                         last_gap  = gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' => ... (%s in tuple)\n'%len(d))
                    break
        f.write(last_gap+")\n")
    elif type(d) == types.DictType:
        if not d:
            f.write(first_gap+"{}\n")
            return
        # recurse on dictionaries
        f.write(first_gap+"{\n")
        keys = d.keys()
        keys.sort()
        key_strings = map(lambda k: ifab(type(k)==types.StringType, k, `k`), keys)
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in map(None, keys, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, d[k],
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if type(k) == types.TupleType:
                            remaining_keys.append(`k`)
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = string.join(remaining_keys,',')
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f, '  %s (%s keys)'%(remaining_keys, len(keys)),0,maxw,0,
                                 gap,gap,'')
                    break
            
            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"}\n")
    elif type(d) == types.InstanceType:
        fields = dir(d)
        
        if not fields:
            f.write(first_gap+"*EmptyClass*\n")
            return
        # recurse on classes
        f.write(first_gap+"*ClassInstance %s\n"%d)
        fields.sort()
        key_strings = map(lambda k: ifab(type(k)==types.StringType, k, `k`), fields)
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in map(None, fields, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, eval('d.'+k),
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if type(k) == type(()):
                            remaining_keys.append(`k`)
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = string.join(remaining_keys,',')
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f,
                                 '  %s (%s keys)'%(remaining_keys, len(keys)),
                                 0,
                                 maxw,
                                 0,
                                 gap,
                                 gap,
                                 '')
                    break
            
            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"*\n")
    elif type(d) == type(""):
        # simply print strings (no quotes)
        if maxw and len(d)>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'..'+d[-final:]+' (%s)\n' % len(d))
        else:
            f.write(first_gap+d+'\n')
    else:
        # string conversion of all other types
        if maxw and len(`d`)>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+`d`[:maxw-final]+'..'+`d`[-final:]+' (%s)\n' % len(`d`))
        else:
            f.write(first_gap+`d`+'\n')

def pp(d,level=-1,maxw=0,maxh=0,parsable=0):
    """ wrapper around pretty_print that prints to stdout"""
    if not parsable: 
        pretty_print(sys.stderr, d, level, maxw, maxh, '', '', '')
    else:
        import pprint
        if maxw: pp2 = pprint.PrettyPrinter(width=maxw, indent=1)#, depth=level
        else: pp2 = pprint.PrettyPrinter(indent=1)#, depth=level
        pp2.pprint(d)
