#!/usr/bin/env python3

import argparse
import json
import os
import re
import sys
import tempfile

# ------------------------------------------------------
# Globals
# ------------------------------------------------------
CSV_FILE_RE = r'^.*\.(csv)$'

# ------------------------------------------------------
# logging
# ------------------------------------------------------
def fatal(msg):
    sys.stderr.write("FATAL - " + msg + "\n")
    sys.exit(1)

def info(msg):
    sys.stderr.write("INFO - " + msg + "\n")

def debug(msg):
    if DEBUG:
        sys.stderr.write("DEBUG - " + msg + "\n")

def warn(msg):
    sys.stderr.write("WARN - " + msg + "\n")

# ------------------------------------------------------
# read_csv_dir
# ------------------------------------------------------
def read_csv_dir(dpath):
    files = []

    for file in sorted(os.listdir(dpath)):
        m = re.match(CSV_FILE_RE, file)
        if m:
            files.append(file)
    info("read files from CSV dir " + dpath + ": " + str(files))
    return files


# ------------------------------------------------------
# add_csv_file
# ------------------------------------------------------
def add_csv_file(js_fh, dpath, file):
    info("add_csv_file " + dpath + " - " + file)
    cpath = os.path.join(dpath, file)

    m = re.match(r'^(\S+)-PAV.*$', file)
    if not m:
        warn("failed to parse sample id from filename " + file)
        m = re.match(r'^(.*).csv$', file)
        
    sample_id = m.group(1)
    info("sample = " + sample_id)
    sample_id = re.sub(r'[\s\-]', '_', sample_id)
    
    # gzip and base64 encode
    enc_cmd = "gzip -c " + cpath + " | uuencode -m " + sample_id + " "
    encoded_data = ''
    
    enc_fh = os.popen(enc_cmd)
    for line in enc_fh:
        if not re.match(r'^(begin-base64.*|====)*$', line):
            encoded_data += line.strip()
        
    # insert encoded data into output JS file
    var_name = "x" + re.sub(r'[\.\-\s]', '_', file) 
    
    js_fh.write("// " + sample_id + "\n")
    js_fh.write("const " + var_name + " = `" + encoded_data + "`;\n\n")

    return { 'sample_id': sample_id, 'var': var_name, 'file': file }
    

# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():
    # input
    parser = argparse.ArgumentParser(description='Compress and uuencode CSV files for single-file distribution.')
    parser.add_argument('--csv_dir', required=True, help='Path to directory that contains the data files to inline.')
    parser.add_argument('--js_output', required=True, help='Path to JavaScript output file.')
    args = parser.parse_args()
    
    files = read_csv_dir(args.csv_dir)

    with open(args.js_output, "w") as js_fh:
        js_fh.write("import { ungzip } from 'pako';\n\n")
        file_info = []
        
        for file in files:
            v = add_csv_file(js_fh, args.csv_dir, file)
            file_info.append(v)

        js_fh.write("const files = [\n")
        for fi in file_info:
            js_fh.write("  { 'sample_id': '" + fi['sample_id'] + "', 'file': '" + fi['file'] + "', 'encoded_data':  " + fi['var'] + " },\n")
        js_fh.write("];\n\n")

        # extractData
        js_fh.write("function extractData(d) {\n"
                    "  const bstring = window.atob(d.encoded_data);\n"
                    "  const barray = Uint8Array.from(bstring, c => c.charCodeAt(0));\n"
                    "  const unc = ungzip(barray)\n"
                    "  var tdec = new TextDecoder('utf-8');\n"
                    "  var data = tdec.decode(unc);\n"
                    "  return data;\n"
                    "}\n\n"
                    )
        
        js_fh.write("files.forEach(f => { f.data = extractData(f) });\n\n")

        js_fh.write("export function getData() { return files; }\n\n")
        
if __name__ == '__main__':
    main()

