import resource
from time import sleep

from pyarrow import Table
from pyarrow import csv
import pyarrow as pa
from pyarrow.parquet import ParquetWriter

import pandas as pd

import click

from memory_profiler import profile

from datetime import timedelta
from timeit import time

import pyarrow as pa
@click.command()
@click.option('--file', prompt='Absolute path to target CSV', help='target csv to extract', required=True)
@click.option('--output_dir', prompt='Absolute path target directory. Use slash to end the path', help='directory to save paritions', required=True)
@click.option('--output_prefix', prompt='The prefix for the outfiles.', help='prefix for partitioned chunks')
@click.option('--column1', prompt='Title of the smiles column (case sensitive)', required=True)
@click.option('--column2', prompt='Title of the name column (case sensitive)', required=True)



@profile
def main(file, output_dir, output_prefix, column1, column2):


    def stopwatch(method):
        def timed(*args, **kw):
            ts = time.perf_counter()
            result = method(*args, **kw)
            te = time.perf_counter()
            duration = timedelta(seconds=te - ts)
            print(f"{method.__name__}: {duration}")
            return result
        return timed


    class InputStreamReader:
        def __init__(self, file_stream):
            self.file_stream = file_stream
            self._stream = None

        def batches(self):
            i = tries = 0
            while True:
                try:
                    batch = self.__next_batch()
                    i += 1
                    yield i, batch
                except StopIteration:
                    break

        def __next_batch(self):
            return self.stream.read_next_batch()


        @property
        def stream(self):
            if not self._stream:
                read_options = pa.csv.ReadOptions(block_size=1048576000)
                parse_options = pa.csv.ParseOptions(delimiter='\t')
                convert_options = pa.csv.ConvertOptions(include_columns=include_columns)
                self._stream = pa.csv.open_csv(
                    self.file_stream, read_options=read_options,
                    parse_options=parse_options,
                    convert_options=convert_options
                
            )
            return self._stream


    @stopwatch        
    def csv_stream_to_parquet_batch_writer(include_columns, \
                     input_file_to_stream, \
                     output_file_stream_directory, \
                     output_file_prefix, \
                     smiles_column_title):
        
        print('Initiating stream.')
        
        input_stream_reader = InputStreamReader(input_file_to_stream)
        
        outfiles_list = []
        
        for i, batch in input_stream_reader.batches():
            print(f'Ingesting batch number {i}')
            df = batch.to_pandas()
            table = pa.Table.from_pandas(df)
            schema = table.schema
            smiles = list(df[smiles_column_title])
            print(f'Writing a total of {len(smiles)} smiles per output file to disk.')
            outfile = f'{output_file_stream_directory}{output_file_prefix}_{i}.parquet'
            ParquetWriter(outfile, schema).write_table(table)
            print(f'Wrote parquet to {outfile}')
            outfiles_list.append(outfile)
            
        return outfile



    include_columns = [column1, column2]

    outfiles = csv_stream_to_parquet_batch_writer(include_columns, \
                         file, \
                         output_dir, \
                         output_prefix, \
                         column1)

    print(outfiles)



if __name__ == '__main__':
    main()





