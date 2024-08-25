import ctypes

# Load the compiled shared library
lib = ctypes.CDLL('/blue/boucher/tyler.pencinger/CUDARHLibrary.so')

# Define the argument types and return type for the processFASTA function
lib.processFASTA.argtypes = [ctypes.c_char_p, ctypes.c_int, ctypes.c_int, ctypes.c_int]
lib.processFASTA.restype = None  # No return value

# Prepare the arguments to pass to the function
file_path = b"/blue/boucher/tyler.pencinger/sequences.fasta"  
window_size = 2  # Example window size
base = 31        # Example base (a prime number)
mod = 1000000007 # Example modulus (a large prime number)

# Call the function from the shared library
lib.processFASTA(file_path, window_size, base, mod)