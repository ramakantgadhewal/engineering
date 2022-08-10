import numpy as np

def get_ideal_number_of_bins_per_sturgis_rule(n=31):
        
    # Sturgis rule for the optimal number of bins in a histogram. Taken from 
    # Bechwith, Mechanical Measurements, 5th edition page 65 footnote *
    # N = 1 + 3.3 log n, n is the total number of points


    N = round(1 + 3.3 * np.log10(n)) # 5.92 or 6

    # https://www.statology.org/sturges-rule/ equation is log2(n) + 1

    N = round(np.log2(n) + 1) # 5.95 or 6

    return N



if __name__ == "__main__":
    
    def main():
        N = get_ideal_number_of_bins_per_sturgis_rule()
        print(f"Ideal number of bins is: {N}")
    main()