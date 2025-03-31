import os
import sys
import csv

def save_rows_as_files(csv_path):
    if not os.path.isfile(csv_path):
        print(f"Error: The file {csv_path} does not exist.")
        return
    
    with open(csv_path, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            if row:
                filename = row[0] + ".csv"
                with open(filename, 'w', newline='') as newfile:
                    csvwriter = csv.writer(newfile)
                    csvwriter.writerow(row)
                print(f"Created file: {filename}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <relative_path_to_csv>")
        return
    
    relative_path = sys.argv[1]
    absolute_path = os.path.abspath(relative_path)
    
    save_rows_as_files(absolute_path)

if __name__ == "__main__":
    main()
