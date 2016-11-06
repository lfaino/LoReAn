#!/usr/bin/env python

''' 
MASTER THESIS PROJECT
Author: Jose A. Espejo
Date: September 2015 - March 2016

Assembling and preparing transcripts
'''
import os 
import subprocess

def trinity(bam_file, wd, max_intron_length, threads):
    '''Calls genome guided trinity on the BAM file to generate 
    assembled transcripts'''
    real_threads = (int(threads))/2
    out_dir = wd+'trinity_out_dir/'
    args = ['Trinity', '--genome_guided_bam', bam_file,
            '--genome_guided_max_intron', max_intron_length, '--max_memory', '10G',
            '--output', out_dir, '--CPU', str(real_threads)]
    out_name = out_dir + 'Trinity-GG.fasta'
    
    if os.path.isfile(out_name): 
        print ('Trinity-GG file existed already: ' + out_name + ' --- skipping\n')
        return out_name
    log_name_err = wd + 'trinity.err.log'
    log_err = open(log_name_err, 'w')
    log_name = wd + 'trinity.log'
    log = open(log_name, 'w')
    try:
        subprocess.check_call(args, stdout = log,stderr = log_err)
    except:
        print 'Trinity did not work properly\n'
        raise NameError('')
    log_err.close()
    log.close()
    return out_name

def seqclean(trinity_file, wd):
    '''the function prepare the Trinity fasta file for PASA '''
    out_dir = '/'.join(trinity_file.split('/')[:-1])+'/'
    t_name = trinity_file.split('/')[-1]
    args = ['seqclean', trinity_file, '-r', out_dir+t_name+'.cln',
            '-o', out_dir+t_name+'.clean']
    out_name = trinity_file + '.clean'
    if os.path.isfile(out_name): 
        print ('Cleaned transcript file existed already: ' + out_name + ' --- skipping\n')
        return out_name
    log_name = wd + 'seqclean.log'
    log = open(log_name, 'w')
    try:
        subprocess.check_call(args, stderr=log, cwd = out_dir)
    except:
        print 'Seqclean did not work properly\n'
        raise NameError('')
    log.close()
    return out_name

def pasa_configuration(pasa_dir, pasa_db):
    '''Creates a PASA configuration file. Database name will be the reference name'''
    conf_file = pasa_dir + 'alignAssembly.config'
    if os.path.isfile(conf_file): 
        print ('PASA configuration file existed already: ' + conf_file + ' --- skipping\n')
        return conf_file
    conf = open(conf_file, 'w')
    lines = ['MYSQLDB='+pasa_db, 
             'validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>',
             'validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>',
             'subcluster_builder.dbi:-m=50']
    for line in lines:
        conf.write(line + '\n')   
    conf.close()
    return conf_file

def pasa_call(pasa_dir, conf_file, pasa_db, reference, transcripts, max_intron_length,threads):
    '''PASA to construct a database of transcripts. It will overwrite any 
    database with the same name -the one of the reference-.'''
    args = ['Launch_PASA_pipeline.pl', '-c', conf_file, '-C', '-r', '-R', '-g',
            reference, '-t', transcripts, '--ALIGNERS', 'gmap', '-I', max_intron_length, '--CPU', str(threads)]
    out_file = pasa_dir + pasa_db + '.pasa_assemblies.gff3'
    #print out_file, os.path.isfile(out_file)
    if os.path.isfile(out_file): 
        print ('PASA output existed already: ' + out_file + ' --- skipping\n')
        return out_file
    log_name = pasa_dir + 'pasa.err.log'
    log_out_name = pasa_dir + 'pasa.out.log'
    log = open(log_name, 'w') 
    out_log = open(log_out_name, 'w')
    try:
        subprocess.check_call(args, stdout = out_log, stderr = log, cwd = pasa_dir)
    except:
        print 'PASA failed'
        raise NameError('')
    log.close()
    out_log.close()
    return out_file

def braker_call(wd, reference, bam_file, species_name, threads, fungus):
    '''Calls braker, may take a while'''
    #perl ~/bin/BRAKER1/braker.pl --cores=3 --workingdir=/home/jose/mapper_testing/braker1_output/ --species=gmap_gff3 
    #--genome=/home/jose/Reference/JR2_Chr8/Verticillium_dahliaejr2.GCA_000400815.2.29.dna.chromosome.8.fa 
    #--bam=/home/jose/mapper_testing/gmap/gmap_Chr8_2Dall.sorted.bam
    if fungus:
        args = ['braker.pl', '--cores='+str(threads), '--useexisting', '--species='+species_name, 
            '--workingdir='+wd, '--genome='+reference, '--fungus' , '--bam='+bam_file]

    else:
        args = ['braker.pl', '--cores='+str(threads), '--useexisting', '--species='+species_name, 
            '--workingdir='+wd, '--genome='+reference, '--bam='+bam_file]
    out_dir = wd+'braker/'
    if os.path.isdir(out_dir): 
        print ('BRAKER1 output existed already: ' + out_dir + ' --- skipping\n')
        return out_dir
    log_name = wd + 'braker.log'
    log = open(log_name, 'w')      
    try:
        braker_ex = subprocess.Popen(args, stderr = log)
        braker_ex.communicate()
    except:
        raise NameError('')
    log.close()   
    return out_dir
def augustus_call(wd, ref, species_name):
    args = ['augustus','--species='+species_name, ref]
    chromo = ref.split('/')[-1]
    wd_augu = wd + '/' + chromo + '.augustus.gff'
    wd_output = wd + '/' + 'augustus.gff3'
    if os.path.isfile(wd_output): 
        print ('Augustus  files exist: ' + wd_output + ' --- skipping\n')
        return  wd
    else:    
        log_name =  wd_augu
        log = open(log_name, 'w')      
        log_name_err = wd + 'augustus.err.log'
        log_e = open(log_name_err, 'w')  
        try:
            subprocess.check_call(args, stderr = log_e, stdout = log, cwd=wd)
        except:
            raise NameError('')
        log.close()   
        log_e.close()
    return wd

def gmes_call(wd, ref, fungus, threads):
    wd_output = wd + '/genemark.gtf.gff3'
    #print wd_output
    if os.path.isfile(wd_output): 
        print ('GeneMarkES  files exist: ' + wd_output + ' --- skipping\n')
        return  wd
    else:
        if fungus:
            args = ['gmes_petap.pl', '--ES', '--fungus' , '--core', threads, '--sequence',ref]
            log_name = wd + 'gmes_petap.gff'
            log = open(log_name, 'w')      
            log_name_err = wd + 'gmes_petaps.err.log'
            log_e = open(log_name_err, 'w')  
            try:
                subprocess.check_call(args, stderr = log_e, stdout = log, cwd=wd)
            except:
                raise NameError('')
            
            log.close()   
            log_e.close()
            
        else:
            args = ['gmes_petap.pl', '--ES', '--core', threads, '--sequence',ref]
            log_name = wd + 'gm_es.gff'
            log = open(log_name, 'w')      
            log_name_err = wd + 'gm_es.err.log'
            log_e = open(log_name_err, 'w')  
            try:
                subprocess.check_call(args, stderr = log_e, stdout = log, cwd=wd)
            except:
                raise NameError('')
            log.close()   
            log_e.close()
            
        return wd
