#!/usr/bin/env python

import os
import re
import shutil
import Image

NORMAL_FRAME_RATE = '25'
SLOW_FRAME_RATE = '1.5'
OUT_FRAME_RATE = '25'

def build_command(commandbuilder, framerate, directory, block_i):
    fn_format_str = "%s/block%d-frame%%d.png" % (directory, block_i)
    commandbuilder.append('ffmpeg')
    commandbuilder.append('-r %s' % framerate)
    commandbuilder.append('-i %s' % fn_format_str)
    commandbuilder.append('-pix_fmt yuv420p')
    commandbuilder.append('-r %s' % OUT_FRAME_RATE)
    commandbuilder.append('%s/inter-block%d.mp4' % (directory, block_i))
    commandbuilder.append('\n')
    
def get_and_sort_filenames(directory, startswith):
    files = os.listdir(directory)
    augfiles = list()
    for f in files:
        if f.startswith(startswith):
            augfiles.append([f, int(re.search(r'\d+', f).group())])
    sortedaugfiles = sorted(augfiles, key=lambda j: j[1])
    
    return sortedaugfiles

def add_to_rename_list(accel_oldname, onetr_oldname, newname, accel_dir, onetr_dir, accel_renamelist, onetr_renamelist):
    accel_renamelist.append((accel_oldname, "%s/%s" % (accel_dir, newname)))
    onetr_renamelist.append((onetr_oldname, "%s/%s" % (onetr_dir, newname)))
    
def concatenate_pngs(accel_old, onetr_old, new_fn):    
    accel_image = Image.open(accel_old)
    onetr_image = Image.open(onetr_old)
    
    width = accel_image.size[0]
    height = accel_image.size[1]
    
    out_image = Image.new("RGB", (2 * width, height))
    
    out_image.paste(accel_image, (0, 0))
    out_image.paste(onetr_image, (width, 0))
    
    out_image.save(new_fn)

def main(accel_dir='accel6-3/', onetr_dir='onetraj3-3/', accel_output_movie='accel6-3.mp4', onetr_output_movie='onetr3-3.mp4', copy=True, combine=True):
    accel_files = get_and_sort_filenames(accel_dir, 'frame')
    onetr_files = get_and_sort_filenames(onetr_dir, 'frame')
    print("Got filenames")
    
    accel_renamelist = list()
    onetr_renamelist = list()
    
    commandbuilder1 = ['\n']
    commandbuilder2 = ['\n']
    normal = True
    block_i = 0
    frame_i = 1
    
    
    
    onetr_fi = 0
    for accel_fi in range(len(accel_files)):
        accel_oldname = accel_files[accel_fi][0]
        onetr_oldname = onetr_files[onetr_fi][0]
        
        
        
        if re.search(r'-cl', accel_oldname) is None:
            # We are in a normal zone
            if normal is None or not normal:
                # We are switching blocks to normal
                
                # Build a command
                build_command(commandbuilder1, framerate, accel_dir, block_i)
                build_command(commandbuilder2, framerate, onetr_dir, block_i)
                
                # Record the change
                normal = True
                block_i += 1
                frame_i = 1
            framerate = NORMAL_FRAME_RATE
            onetr_fi += 1
            
        else:
            # We are in an abnormal zone
            if normal is None or normal:
                # We are switching blocks to slow
                
                # Build a command
                build_command(commandbuilder1, framerate, accel_dir, block_i)
                build_command(commandbuilder2, framerate, onetr_dir, block_i)        
                
                # Record the change
                normal = False
                block_i += 1
                frame_i = 1
                
                if re.search(r'-cl\.', accel_oldname) is not None:
                    # For some reason, it drops the first frame, so we'll make two
                    newname = "block%d-frame%d.png" % (block_i, frame_i)
                    add_to_rename_list(accel_oldname, onetr_oldname, newname, accel_dir, onetr_dir, accel_renamelist, onetr_renamelist)
                    frame_i += 1
            
            framerate = SLOW_FRAME_RATE

        # Get ready to rename the files                
        newname = "block%d-frame%d.png" % (block_i, frame_i)
        add_to_rename_list(accel_oldname, onetr_oldname, newname, accel_dir, onetr_dir, accel_renamelist, onetr_renamelist)
        
        
        frame_i += 1

    print("Generated rename list")
    
    # Now add one command to string all the intermediates together
    commandbuilder1.append('ffmpeg')
    commandbuilder2.append('ffmpeg')
    
    for i in range(block_i):
        fn_str = '%s/inter-block%d.mp4' % (accel_dir, i)
        commandbuilder1.append('-i %s' % fn_str)
        
        fn_str = '%s/inter-block%d.mp4' % (onetr_dir, i)
        commandbuilder2.append('-i %s' % fn_str)
    
    commandbuilder1.append('-filter_complex concat=n=%d' % (block_i))
    commandbuilder1.append('%s' % accel_output_movie)
    commandbuilder1.append('\n')
    
    commandbuilder2.append('-filter_complex concat=n=%d' % (block_i))
    commandbuilder2.append('%s' % onetr_output_movie)
    commandbuilder2.append('\n')
    
    commandbuilder1.extend(commandbuilder2)
    
    command = ' '.join(commandbuilder1)
    
    
    with open('makemovie.sh', 'w') as sh:
        sh.write(command)
        
    print("Wrote command")
    

    
    
    assert len(accel_renamelist) == len(onetr_renamelist)
    
    if combine:
        for t in range(len(accel_renamelist)):
            old_accel = '%s/%s' % (accel_dir, accel_renamelist[t][0])
            old_onetr = '%s/%s' % (onetr_dir, onetr_renamelist[t][0])
            new = '%s' % (onetr_renamelist[t][1])
            concatenate_pngs(old_accel, old_onetr, new)
            if copy:
                # Copy just the accelerated one as well
                new = accel_renamelist[t][1]
                shutil.copy(old_accel, new)
    else:
        print("Copying/renaming files 1")
        for t in accel_renamelist:
            old = '%s/%s' % (accel_dir, t[0])
            new = '%s' % (t[1])
            if copy:
                shutil.copy(old, new)
                
        print("Copying/renaming files 2")
                
        for t in onetr_renamelist:
            old = '%s/%s' % (onetr_dir, t[0])
            new = '%s' % (t[1])
            if copy:
                shutil.copy(old, new)
            
    print("Done")
    
