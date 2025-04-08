from python.segmentation.train import ManagedList

class TestManagedList():

    # Initialize ManagedList with values and verify they are stored correctly
    def test_initialize_with_values(self):
        # Given
        initial_values = [1, 2, 3, 4, 5]
        max_size = 10
    
        # When
        managed_list = ManagedList(initial_values, max_size)
    
        # Then
        assert managed_list.get_list() == initial_values
        assert len(managed_list.values) == len(initial_values)
        assert managed_list.values.maxlen == max_size

    # Initialize with empty list and verify behavior
    def test_initialize_with_empty_list(self):
        # Given
        empty_list = []
    
        # When
        managed_list = ManagedList(empty_list)
    
        # Then
        assert managed_list.get_list() == []
        assert len(managed_list.values) == 0
        assert managed_list.mean() is None
        assert managed_list.median() is None

    # Add multiple values in sequence and verify order is maintained
    def test_add_multiple_values_maintains_order(self):
        # Given
        initial_values = [1, 2, 3]
        new_values = [4, 5, 6]
        expected_values = [1, 2, 3, 4, 5, 6]
        managed_list = ManagedList(initial_values)
    
        # When
        for value in new_values:
            managed_list.add(value)
    
        # Then
        assert managed_list.get_list() == expected_values

    # Add more values than max size
    def test_add_more_than_max_size(self):
        # Given
        initial_values = [1, 2, 3, 4, 5]
        max_size = 5
        managed_list = ManagedList(initial_values, max_size)
    
        # When
        managed_list.add(6)
        managed_list.add(7)
    
        # Then
        assert managed_list.get_list() == [3, 4, 5, 6, 7]
        assert len(managed_list.values) == max_size