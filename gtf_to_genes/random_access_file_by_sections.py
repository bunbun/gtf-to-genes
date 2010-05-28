#!/usr/bin/env python
"""
    Divide file into sections and write file position directory to each section
    This allows entire sections to be skipped
"""
import struct


#_____________________________________________________________________________________
#
#   read_directory_of_sections
#_____________________________________________________________________________________
def read_directory_of_sections (stream):
    """
    Random access to different sections
    """
    cnt_sections = struct.unpack("q", stream.read(8))[0]

    file_pos_by_section = dict()
    for i in range(cnt_sections):
        section_name_len = struct.unpack("q", stream.read(8))[0]
        section_name     = stream.read(section_name_len)
        section_file_pos = struct.unpack("q", stream.read(8))[0]
        file_pos_by_section[section_name] = section_file_pos

    return file_pos_by_section

#_____________________________________________________________________________________
#
#   write_directory_of_sections
#_____________________________________________________________________________________
def write_directory_of_sections (stream, sections):
    """
    Random access to different sections
    Create directory for each of the different random access sections
    The file positions will be filled in once we have written the actual sections
    """
    # save positions in section directory which will be filled later
    section_name_directory_entry_pos = dict()


    # number of sections
    stream.write(struct.pack("q", len(sections)))


    for section_name in sections:
        # name of section_name
        stream.write(struct.pack("q", len(section_name)))
        stream.write(section_name)

        # where the section starts: currently only placeholder
        section_name_directory_entry_pos[section_name] = stream.tell()
        stream.write(struct.pack("q", 0))

    return section_name_directory_entry_pos

#_____________________________________________________________________________________
#
#   fill_directory_of_sections
#_____________________________________________________________________________________
def fill_directory_of_sections (stream, section_name_directory_entry_pos, file_pos_by_section):
    """
    Random access to different sections
    Fill positions in directory for each of the different random access sections
    """
    for section_name, section_file_pos in file_pos_by_section.iteritems():
        directory_entry_pos = section_name_directory_entry_pos[section_name]
        stream.seek(directory_entry_pos)
        stream.write(struct.pack("q", section_file_pos))








