import bisect
import numpy as np

class OutlierList:
    def __init__(self, values : list=[], max_size : int=100):
        self.values = []  # Sorted list to maintain the values
        self.max_size = max_size
        
        for v in values:
            self.append(v)

    def _iqr(self):
        """Helper method to calculate IQR (Interquartile Range)."""
        if len(self.values) < 4:
            return None, None, None  # Not enough data for IQR
        q1 = np.percentile(self.values, 25)
        q3 = np.percentile(self.values, 75)
        iqr = q3 - q1
        return q1, q3, iqr

    def append(self, new_value):
        """Add a new value, replacing an outlier if needed."""
        if len(self.values) < self.max_size:
            # List has fewer than max_size, just insert in the right position
            bisect.insort(self.values, new_value)
        else:
            # Calculate IQR to identify outliers
            q1, q3, iqr = self._iqr()
            if iqr is not None:
                lower_bound = q1 - 1.5 * iqr
                upper_bound = q3 + 1.5 * iqr

                # If the new value is within bounds
                if lower_bound <= new_value <= upper_bound:
                    # Replace the most extreme outlier or the oldest value
                    outliers = [v for v in self.values if v < lower_bound or v > upper_bound]
                    if outliers:
                        # Remove the most extreme outlier
                        if abs(outliers[0] - new_value) > abs(outliers[-1] - new_value):
                            self.values.remove(outliers[0])
                        else:
                            self.values.remove(outliers[-1])
                    else:
                        # No outliers: replace the oldest value (first in the sorted list)
                        self.values.pop(0)
                    bisect.insort(self.values, new_value)

    def mean(self):
        return np.mean(self.values)
    
    def median(self):
        return np.median(self.values)
    
    def std(self):
        return np.std(self.values)

    def get_values(self):
        """Return the current list of values."""
        return self.values