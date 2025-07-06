
"""
author: Hadi Vareno
e-mail: mohammad.noori.vareno@uni-jena.de
github: https://github.com/TheVareno
"""
 
from read5.Reader import read # type: ignore 
# import ont_fast5_api   # type: ignore 
from ont_fast5_api.conversion_tools.fast5_subset import Fast5Filter   # type: ignore 
import argparse
# TODO use hampelFilter from FileIO.py
from hampel import hampel # type: ignore
import subprocess as sp
import os  
import multiprocessing as mp 
import queue



def setup_working_directory(working_dir=None)-> None:
    if working_dir is None:
        working_dir = os.getcwd()
    if not os.path.exists(working_dir):
        raise FileNotFoundError(f"The working directory: '{working_dir}' does not exist.")
    
    os.chdir(working_dir)


    
def find_polya(task_queue: mp.Queue, result_queue: mp.Queue, read_object: str): 
    
    while not task_queue.empty(): 
        try:
            read_id = task_queue.get_nowait() 
            z_normalized_signal_values = read_object.getZNormSignal(read_id, mode='mean')
            filter_object = hampel(z_normalized_signal_values, window_size=5, n_sigma=6.0)
            filtered_signal_values = filter_object.filtered_data

            if len(filtered_signal_values) == 0:
                print(f"the array of signal values empty for read id : {read_id}")
                
            polyA_app_call = './polyA'  
            
            sig_vals_str = ','.join(map(str, filtered_signal_values))
            
            process = sp.Popen(polyA_app_call, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
            
            if not sig_vals_str: 
                print(f"Empty signal values for read {read_id}")
            
            process.stdin.write(f"{sig_vals_str}\n")
            process.stdin.flush()
            stdout, stderr = process.communicate()
            rc = process.returncode # returns int 
            
            if rc == 0:  
                borders = stdout.strip()
                result_queue.put((read_id, borders))
            else: 
                pass
            
            if stderr: 
                print(f"Error for {read_id}: {stderr}")
                continue

        except queue.Empty:
            break



def split_segment_input(input_read_data: str, output_path: str, summary_file_path: str):

    raw_signal_path = input_read_data  

    if not os.path.exists(output_path): # dir exist. check 
        os.makedirs(output_path) 
    
    save_file = os.path.join(output_path, 'region_borders.csv')  
    with open(save_file, 'w') as f: # file exist. check 
        f.write("Read ID,Poly(A) end,Adapter end,Leader end,Start end\n")
    
    splitter = Fast5Filter(input_folder=raw_signal_path, 
                output_folder=output_path, 
                read_list_file=summary_file_path,
                filename_base="subset",
                batch_size=500, 
                threads=1,
                recursive=False,
                file_list_file=None,
                follow_symlinks=False,
                target_compression=None)

    splitter.run_batch() # raw signal splitting into 8 files of 500 reads 

    for file in os.listdir(output_path):
        filename = os.fsencode(file)
        
        if filename.endswith(b".fast5") or filename.endswith(b".pod5") or filename.endswith(b".slow5"): 
            
            file = os.path.join(output_path, file)
            read_object = read(file) # needs file path ends with .fast5 / .pod5 / .slow5
            
            all_read_ids = read_object.getReads() # 500 each time
            
            task_queue = mp.Queue()
            result_queue = mp.Queue()
    
            for read_id in all_read_ids:
                task_queue.put(read_id)
    
            number_of_processes = os.cpu_count()
        
            processes = [ mp.Process(target=find_polya, args=(task_queue, result_queue, read_object)) for _ in range(number_of_processes) ]
    
            for proc in processes:
                proc.start()
        
            for proc in processes:
                proc.join()
    
            while not result_queue.empty():
                read_id, borders = result_queue.get()
                with open (save_file, 'a') as f: 
                    f.write(f"{read_id},{borders}\n")
        else: 
            continue


def main()-> None: 
    # TODO use same parameters/namings as dynamont
    parser = argparse.ArgumentParser(description="Process and Save output file.")
    parser.add_argument('-i', "--fast5_path", type=str, required=True, help="Path to input ONT read data in FAST5, POD5, or SLOW5 format.")
    parser.add_argument('-o', "--output_dir", type=str, required=True, help="Directory to save output files.")
    parser.add_argument('-s', "--summary_file", type=str, required=True, help="Path to the sequence summary file.")
    parser.add_argument('-w', "--working_dir", type=str, default=os.getcwd(), help="Working directory for the program (default is the current directory).")

    args = parser.parse_args()

    try:
        setup_working_directory(working_dir=args.working_dir)
    
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    split_segment_input(input_read_data=args.fast5_path, 
                        output_path=args.output_dir, 
                        summary_file_path=args.summary_file)

    

if __name__ == '__main__' : 
    main()






