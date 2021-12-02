import calendar
import time

sort_by_index: int = 4
header_line = "source_id\tra_ep2000\tdec_ep2000\tphot_g_mean_mag\tdistance\n"


def generate_output_file(list_of_N_stars: list,
                         output_file: str,
                         separator: str,
                         end_of_line: str):
    """ Open the file for writing the data of N stars to the output file contained in the field of view. """

    with open(output_file, 'w') as fw:
        fw.write(header_line)
        for i in range(len(list_of_N_stars)):
            star_components: str = ''
            for j in range(len(list_of_N_stars[i])):
                if j < len(list_of_N_stars[i]) - 1:
                    star_components += str(list_of_N_stars[i][j]) + separator
                else:
                    star_components += str(list_of_N_stars[i][j]) + end_of_line
            fw.write(star_components)


def sort_by_distance(list_of_stars_inside_the_fov: list):
    """ Function to sort the list by distance of the stars from a given point.

    :param list_of_stars_inside_the_fov: Stars contained inside the field of view
    :type list_of_stars_inside_the_fov: list
    """

    if len(list_of_stars_inside_the_fov) > 1:
        mid = len(list_of_stars_inside_the_fov) // 2
        left_half = list_of_stars_inside_the_fov[:mid]
        right_half = list_of_stars_inside_the_fov[mid:]

        sort_by_distance(left_half)
        sort_by_distance(right_half)

        i = j = k = 0
        while i < len(left_half) and j < len(right_half):
            if left_half[i][sort_by_index] < right_half[j][sort_by_index]:
                list_of_stars_inside_the_fov[k] = left_half[i]
                i += 1
            else:
                list_of_stars_inside_the_fov[k] = right_half[j]
                j += 1
            k += 1

        while i < len(left_half):
            list_of_stars_inside_the_fov[k] = left_half[i]
            i += 1
            k += 1

        while j < len(right_half):
            list_of_stars_inside_the_fov[k] = right_half[j]
            j += 1
            k += 1


def check_if_star_is_in_fov(dict_of_stars: dict,
                            header_name_ra: str,
                            header_name_dec: str,
                            left_border: float,
                            right_border: float,
                            bottom_border: float,
                            top_border: float):

    ra_of_star = float(dict_of_stars[header_name_ra])
    dec_of_star = float(dict_of_stars[header_name_dec])

    if left_border < ra_of_star < right_border \
            and bottom_border < dec_of_star < top_border:
        return True
    else:
        return False


def calculate_distance(dict_of_stars: dict,
                       header_name_ra: str,
                       header_name_dec: str,
                       ra_coordinate: float,
                       dec_coordinate: float):

    ra_of_star = float(dict_of_stars[header_name_ra])
    dec_of_star = float(dict_of_stars[header_name_dec])
    return ((ra_of_star - ra_coordinate) ** 2 + (dec_of_star - dec_coordinate) ** 2) ** 0.5


def find_n_brightest_stars(file_of_data: str,
                           output_file: str,
                           ra_coordinate: float,
                           dec_coordinate: float,
                           fov_horizontal_length: float,
                           fov_vertical_length: float,
                           number_of_brightest_stars: int,
                           separator: str,
                           end_of_line: str,
                           header_name_ra: str,
                           header_name_dec: str,
                           header_name_source_id: str,
                           header_name_magnitude: str):
    """ Function to find the brightest N stars, sort it by their distance from a given point and write their data
    into the output file.

    :param file_of_data: File to read the data
    :type file_of_data: str
    :param output_file: File for writing the output of the program
    :type output_file: str
    :param ra_coordinate: First coordinate of the given point
    :type ra_coordinate: float
    :param dec_coordinate: Second coordinate of the given point
    :type dec_coordinate: float
    :param fov_horizontal_length: Horizontal field of view length
    :type fov_horizontal_length: float
    :param fov_vertical_length: Vertical field of view length
    :type fov_vertical_length: float
    :param number_of_brightest_stars: Number of stars that should be extracted
    :type number_of_brightest_stars: int
    :param separator: by which the words shall be separated
    :type: str
    :param end_of_line: New line to add in the end of line
    :type end_of_line: str
    :param header_name_ra: Star's first parameter
    :type header_name_ra: str
    :param header_name_dec: Star's second parameter
    :type header_name_dec: str
    :param header_name_source_id: Star's third parameter
    :type header_name_source_id: str
    :param header_name_magnitude: Star's fourth
    :type header_name_magnitude: str
    """

    left_border = ra_coordinate - fov_horizontal_length / 2
    right_border = ra_coordinate + fov_horizontal_length / 2
    bottom_border = dec_coordinate - fov_vertical_length / 2
    top_border\
        = dec_coordinate + fov_vertical_length / 2
    temp_dict = {}
    list_of_stars_in_the_fov = []
    count = 0

    with open(file_of_data, 'r') as fr:
        for line in fr.readlines()[1:]:
            """ Open the file and start reading, check whether the star is in the field of view and if yes, 
            calculate the star's distance from the center of the field of view (given point) and add first 
            N (number_of_brightest_stars) stars to the list (list_of_stars_in_the_fov) in ascending order
            by their magnitude. """

            count += 1
            if count == 1:
                first_line = line.split()
                continue
            line = [line.split(separator) for _ in line.splitlines()][0]
            for i in range(len(first_line)):
                temp_dict[first_line[i]] = line[i]
            star_is_in_fov = check_if_star_is_in_fov(temp_dict,
                                                     header_name_ra,
                                                     header_name_dec,
                                                     left_border,
                                                     right_border,
                                                     bottom_border,
                                                     top_border)
            if star_is_in_fov:
                distance = calculate_distance(temp_dict,
                                              header_name_ra,
                                              header_name_dec,
                                              ra_coordinate,
                                              dec_coordinate)
                temp_list = list([temp_dict[header_name_source_id],
                                  float(temp_dict[header_name_ra]),
                                  float(temp_dict[header_name_dec]),
                                  float(temp_dict[header_name_magnitude]),
                                  distance])
                if len(list_of_stars_in_the_fov) == 0:
                    list_of_stars_in_the_fov.append(temp_list)
                    continue
                for index in range(len(list_of_stars_in_the_fov)):
                    if temp_list[3] < list_of_stars_in_the_fov[index][3]:
                        list_of_stars_in_the_fov.insert(index, temp_list)
                        if len(list_of_stars_in_the_fov) > number_of_brightest_stars:
                            list_of_stars_in_the_fov.pop()
                        break
    sort_by_distance(list_of_stars_in_the_fov)
    generate_output_file(list_of_stars_in_the_fov, output_file, separator, end_of_line)


def check_float(potential_float: str):
    """ Function to check whether the inputted string is a float number and throw an exception if not. """

    try:
        float_number = float(potential_float)
        return float_number
    except ValueError:
        return potential_float


def check_int(potential_int: str):
    """ Function to check whether the inputted string is an integer number and throw an exception if not. """

    try:
        int_number = int(potential_int)
        return int_number
    except ValueError:
        return potential_int


def main():
    """ Driver code. """

    input_file: str = 'cleaned_stars.tsv'
    output_file: str = str(calendar.timegm(time.gmtime())) + '.csv'
    header_name_ra: str = 'ra_ep2000'
    header_name_dec: str = 'dec_ep2000'
    header_name_source_id: str = 'source_id'
    header_name_magnitude: str = 'phot_g_mean_mag'
    separator: str = '\t'
    end_of_line: str = '\n'
    ra_value = check_float(input('Input ra_value: '))
    while type(ra_value) == str:
        print('ValueError!')
        ra_value = check_float(input('Input ra_value: '))
    dec_value = check_float(input('Input dec_value: '))
    while type(dec_value) == str:
        print('ValueError!')
        dec_value = check_float(input('Input dec_value: '))
    horizontal_fov_value = check_float(input('Input horizontal field of view: '))
    while type(horizontal_fov_value) == str:
        print('ValueError!')
        horizontal_fov_value = check_float(input('Input horizontal_fov_value: '))
    vertical_fov_value = check_float(input('Input vertical field of view: '))
    while type(vertical_fov_value) == str:
        print('ValueError!')
        vertical_fov_value = check_float(input('Input vertical_fov_value: '))
    number_of_stars = check_int(input('Input number of stars: '))
    while type(number_of_stars) != int:
        print('ValueError!')
        number_of_stars = check_int(input('Input number_of_stars: '))
    find_n_brightest_stars(input_file,
                           output_file,
                           ra_value,
                           dec_value,
                           horizontal_fov_value,
                           vertical_fov_value,
                           number_of_stars,
                           separator,
                           end_of_line,
                           header_name_ra,
                           header_name_dec,
                           header_name_source_id,
                           header_name_magnitude)


if __name__ != '__main__':
    pass
else:
    main()
