/* ********************************** TREADHDF5 ********************************** */
/*////////////////////////////////////////////////////////////////////////////*/
// PROJECT :        PATHFINDER
// TITLE :          TReadHDF5
// DESCRIPTION :    Utilities to read HDF5 files
// AUTHOR(S) :      Michel-Andrès Breton (michel-andres.breton@obspm.fr)
// CONTRIBUTIONS :  [Michel-Andrès Breton (2015-2021)]
// LICENSE :        CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/
/// \file           TReadHDF5.h
/// \brief          Utilities to read HDF5 files
/// \author         Michel-Andrès Breton (michel-andres.breton@obspm.fr)
/// \date           2015-2021
/// \copyright      CECILL-B License
/*////////////////////////////////////////////////////////////////////////////*/

#ifndef TREADHDF5_H_INCLUDED
#define TREADHDF5_H_INCLUDED

// Include C++
#include <algorithm>
#include <array>
#include <dirent.h>
#include <errno.h>
#include <iostream>
#include <random>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

#ifdef GCCBELOW7
#include <experimental/algorithm>
#endif

#include "hdf5.h"
#include "magrathea/constants.h"

class TReadHDF5 {

    // Methods
public:
    // Native conversion
    static hid_t h5t_native(const int N);
    static hid_t h5t_native(const unsigned int N);
    static hid_t h5t_native(const float N);
    static hid_t h5t_native(const double N);
    static hid_t h5t_native(const unsigned long int N);
    template <typename T>
    static hid_t h5t_native(void);

    // Some informations
    static void pos_from_index_fullsky(const std::string &str, const std::vector<int> &nctab, const std::vector<float> &dimtab, const float &cubesize, std::array<double, 3> &point111);
    static void pos_from_index_narrow(const std::string &str, const std::vector<int> &nctab, const std::vector<float> &dimtab, const float &cubesize, const double &thetay, const double &thetaz, std::array<double, 3> &point111);

    void displayAllInfos(const std::string &fileName) const;
    // Get data
    template <typename Type>
    static void get_data_from_dataset(const hid_t &gid, const std::string &output_name, std::vector<Type> &output);
    template <typename Type, typename... String, typename... Vector>
    static void get_data_from_dataset(const hid_t &gid, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs);
    // Count
    template <typename Integer>
    static void cellsPerLevels(const std::string &fileName, std::vector<Integer> &count);
    template <class Parameter, typename Integer, class Conic>
    static void cellsAndCubesPerLevels(const Parameter &parameters, const std::string &fileName, std::vector<Integer> &count, std::vector<std::vector<std::string>> &cubeNumber, const double &thetay, const double &thetaz, const Conic &conic);

    // Attributes at file level
    template <typename Type>
    static void getAttribute(const std::string &fileName, const std::string &attributeName, Type &attributeValue);
    template <typename Type>
    static void getAttribute(const std::string &fileName, const std::string &attributeName, std::vector<Type> &attributeValue);

    // Attributes on group
    template <typename Type>
    static void getAttribute(const std::string &fileName, const std::string &attributeName1, const std::string &attributeName2, Type &attributeValue);
    template <typename Type>
    static void getAttribute(const std::string &fileName, const std::string &attributeName1, const std::string &attributeName2, std::vector<Type> &attributeValue);

    // Read and fill gravity cell vectors
    // Full sample
    template <typename Type>
    static void fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::string &output_name, std::vector<Type> &output);
    template <typename Type, typename... String, typename... Vector>
    static void fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs);
    // Selected sample
    template <typename Type>
    static void fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::vector<std::string> &cubeNumber, const std::string &output_name, std::vector<Type> &output);
    template <typename Type, typename... String, typename... Vector>
    static void fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::vector<std::string> &cubeNumber, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs);
    // Read and fill particle vectors
    // Full sample
    template <typename Type>
    static void fillVectors_part(const std::string &fileName, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output);
    template <typename Type, typename... String, typename... Vector>
    static void fillVectors_part(const std::string &fileName, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs);
    // Random sample
    template <typename Type>
    static void fillVectors_part(const double &fraction, const std::string &fileName, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output);
    template <typename Type1, typename Type2>
    static void fillVectors_part(const double &fraction, const std::string &fileName, const std::string &fileSide, const std::string &output_name1, std::vector<Type1> &output1, const std::string &output_name2, std::vector<Type2> &output2);
    template <typename Type1, typename Type2, typename Type3>
    static void fillVectors_part(const double &fraction, const std::string &fileName, const std::string &fileSide, const std::string &output_name1, std::vector<Type1> &output1, const std::string &output_name2, std::vector<Type2> &output2, const std::string &output_name3, std::vector<Type3> &output3);
    // Selected sample
    template <class Parameter, typename Type, class Conic>
    static void fillVectors_part(const Parameter &parameters, const std::string &fileName, const double &thetay, const double &thetaz, const Conic &conic, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output);
    template <class Parameter, typename Type, typename... String, typename... Vector, class Conic>
    static void fillVectors_part(const Parameter &parameters, const std::string &fileName, const double &thetay, const double &thetaz, const Conic &conic, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs);
};

/// \brief          Converts Type to the corresponding hid_t.
/// \details        The HDF5 C routines use a hid_t type in order
///                 to get the data of a specified type.
/// \param[in]      N (only the type matters, not the value)
/// \return         hid_t H5T_NATIVE of the intput type
hid_t TReadHDF5::h5t_native(const int N) {
    return H5T_NATIVE_INT;
}

/// \brief          Converts Type to the corresponding hid_t.
/// \details        The HDF5 C routines use a hid_t type in order
///                 to get the data of a specified type.
/// \param[in]      N (only the type matters, not the value)
/// \return         hid_t H5T_NATIVE of the intput type
hid_t TReadHDF5::h5t_native(const unsigned int N) {
    return H5T_NATIVE_UINT;
}

/// \brief          Converts Type to the corresponding hid_t.
/// \details        The HDF5 C routines use a hid_t type in order
///                 to get the data of a specified type.
/// \param[in]      N (only the type matters, not the value)
/// \return         hid_t H5T_NATIVE of the intput type
hid_t TReadHDF5::h5t_native(const float N) {
    return H5T_NATIVE_FLOAT;
}

/// \brief          Converts Type to the corresponding hid_t.
/// \details        The HDF5 C routines use a hid_t type in order
///                 to get the data of a specified type.
/// \param[in]      N (only the type matters, not the value)
/// \return         hid_t H5T_NATIVE of the intput type
hid_t TReadHDF5::h5t_native(const double N) {
    return H5T_NATIVE_DOUBLE;
}

/// \brief          Converts Type to the corresponding hid_t.
/// \details        The HDF5 C routines use a hid_t type in order
///                 to get the data of a specified type.
/// \param[in]      N (only the type matters, not the value)
/// \return         hid_t H5T_NATIVE of the intput type
hid_t TReadHDF5::h5t_native(const unsigned long int N) {
    return H5T_NATIVE_ULONG;
}

/// \brief          Converts Type to the corresponding hid_t.
/// \details        The HDF5 C routines use a hid_t type in order
///                 to get the data of a specified type.
/// \tparam         T Type
/// \return         hid_t H5T_NATIVE of the intput type
template <class T>
hid_t TReadHDF5::h5t_native() {
    return h5t_native(T());
}

/// \brief          Get position from index.
/// \details        Get position from indexes in Raygal simulation fullsky data cubes
/// \param[in]      str Cube name in HDF5 file.
/// \param[in]      nctab Number of cubes in data
/// \param[in]      dimtab Dimensions of data.
/// \param[in]      cubesize Cube size.
/// \param[in,out]  point111 Cube center position.
void TReadHDF5::pos_from_index_fullsky(const std::string &str, const std::vector<int> &nctab, const std::vector<float> &dimtab, const float &cubesize, std::array<double, 3> &point111) {

    std::array<double, 3> point;
    int index(std::stoi(str.substr(4, str.size())));
    const double ihx = int(dimtab[0] / cubesize + 1) * cubesize;
    const double ihy(ihx), ihz(ihx);
    const int iz1 = index % (nctab[0] * nctab[1]) != 0 ? index / (nctab[0] * nctab[1]) + 1 : index / (nctab[0] * nctab[1]);
    index = index % (nctab[0] * nctab[1]);
    const int iy1 = index % (nctab[0] * nctab[1]) == 0 ? (dimtab[0] + ihy) / cubesize + 1 : index % nctab[0] != 0 ? index / nctab[0] + 1
                                                                                                                  : index / nctab[0];
    index = index % nctab[0];
    const int ix1 = index == 0 ? (dimtab[0] + ihx) / cubesize + 1 : index;
    point[0] = (ix1 - 1) * cubesize - ihx;
    point[1] = (iy1 - 1) * cubesize - ihy;
    point[2] = (iz1 - 1) * cubesize - ihz;
    point111[0] = point[0] + cubesize * 0.5;
    point111[1] = point[1] + cubesize * 0.5;
    point111[2] = point[2] + cubesize * 0.5;
}

/// \brief          Get position from index.
/// \details        Get position from indexes in Raygal simulation narrow data cubes
/// \param[in]      str Cube name in HDF5 file.
/// \param[in]      nctab Number of cubes in every dimension.
/// \param[in]      dimtab Dimensions in Ramses units of the
///		    parallelepiped around the cone
/// \param[in]      cubesize Cube size.
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
/// \param[in,out]  point111 Cube center position.
void TReadHDF5::pos_from_index_narrow(const std::string &str, const std::vector<int> &nctab, const std::vector<float> &dimtab, const float &cubesize, const double &thetay, const double &thetaz, std::array<double, 3> &point111) {

    std::array<double, 3> point;
    int index(std::stoi(str.substr(4, str.size())));
    const double hy = dimtab[0] * std::tan(thetay);
    const double hz = dimtab[0] * std::tan(thetaz);
    const double ihy = int((hy / cubesize) + 1) * cubesize;
    const double ihz = int((hz / cubesize) + 1) * cubesize;
    const int iz1 = index % (nctab[0] * nctab[1]) != 0 ? index / (nctab[0] * nctab[1]) + 1 : index / (nctab[0] * nctab[1]);
    index = index % (nctab[0] * nctab[1]);
    const int iy1 = index % (nctab[0] * nctab[1]) == 0 ? (dimtab[0] + ihy) / cubesize + 1 : index % nctab[0] != 0 ? index / nctab[0] + 1
                                                                                                                  : index / nctab[0];
    index = index % nctab[0];
    const int ix1 = index == 0 ? dimtab[0] / cubesize + 1 : index;
    point[0] = (ix1 - 1) * cubesize;
    point[1] = (iy1 - 1) * cubesize - ihy;
    point[2] = (iz1 - 1) * cubesize - ihz;
    point111[0] = point[0] + cubesize * 0.5;
    point111[1] = point[1] + cubesize * 0.5;
    point111[2] = point[2] + cubesize * 0.5;
}

/// \brief          Prints informations.
/// \details        Prints the following informations on the
///		    HDF5 file :
///		    - File name
///		    - Group names
///		    - Number of subgroups
/// \param[in]      fileName File list.
void TReadHDF5::displayAllInfos(const std::string &fileName) const {
    std::cout << "Fichier : " << fileName << std::endl;
    unsigned int const MAX_NAME = 1024;
    char memb_name[MAX_NAME];
    H5G_info_t grpinfo, lvlinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, "/data/", H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    for (unsigned int i = 0; i < grpinfo.nlinks; i++) {
        H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, memb_name, (size_t)MAX_NAME, H5P_DEFAULT);
        std::string levelName(memb_name);
        if (levelName.substr(0, 5) == "level") {
            hid_t lvlid = H5Gopen(gid, memb_name, H5P_DEFAULT);
            H5Gget_info(lvlid, &lvlinfo);
            std::cout << levelName << " (" << lvlinfo.nlinks << " cubes)" << std::endl;
        } //  strcmp
        else {
            std::cout << memb_name << std::endl;
        }
    } //  i
    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from dataset.
/// \details        Get data from dataset
/// \tparam         Type Type type.
/// \param[in]      gid Group HDF5 id.
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
template <typename Type>
void TReadHDF5::get_data_from_dataset(const hid_t &gid, const std::string &output_name, std::vector<Type> &output) {

    hid_t dsid = H5Dopen(gid, output_name.c_str(), H5P_DEFAULT);
    hsize_t alloc = H5Dget_storage_size(dsid);
    unsigned int nmax = alloc / (sizeof(Type));
    Type *dset_data;
    dset_data = (Type *)malloc(sizeof(Type) * nmax);
    H5Dread(dsid, h5t_native<Type>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
    output.insert(output.end(), &dset_data[0], &dset_data[nmax]);
    free(dset_data);
    H5Dclose(dsid);
}

/// \brief          Get data from dataset.
/// \details        Get data from dataset
/// \tparam         Type Type type.
/// \param[in]      gid Group HDF5 id.
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
/// \param[in]	    output_names Variadic Data names.
/// \param[in,out]  outputs Variadic Data vectors.
template <typename Type, typename... String, typename... Vector>
void TReadHDF5::get_data_from_dataset(const hid_t &gid, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs) {

    TReadHDF5::get_data_from_dataset(gid, output_name, output);
    TReadHDF5::get_data_from_dataset(gid, output_names..., outputs...);
}

/// \brief          Number of cells.
/// \details        Number of cells per level
/// \tparam         Integer Count type.
/// \param[in]      fileName File name.
/// \param[in,out]  count Number of cells per level.
template <typename Integer>
void TReadHDF5::cellsPerLevels(const std::string &fileName, std::vector<Integer> &count) {

    unsigned int const MAX_NAME = 1024;
    char memb_name[MAX_NAME];
    H5G_info_t grpinfo, lvlinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, "/data/", H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    // Loop over groups in '/data'
    for (unsigned int i = 0; i < grpinfo.nlinks; i++) {
        H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, memb_name, (size_t)MAX_NAME, H5P_DEFAULT);
        std::string levelName(memb_name);
        // Enter in groups which start with 'level'
        if (levelName.substr(0, 5) == "level") {
            hid_t lvlid = H5Gopen(gid, memb_name, H5P_DEFAULT);
            H5Gget_info(lvlid, &lvlinfo);
            if (lvlinfo.nlinks != 0) {
                count.push_back(0);
                // Get number of cells per level
                hid_t dsid = H5Dopen(lvlid, "ncell_level", H5P_DEFAULT);
                hsize_t alloc = H5Dget_storage_size(dsid);
                unsigned int nmax = alloc / (sizeof(Integer));
                Integer *dset_data;
                dset_data = (Integer *)malloc(sizeof(Integer) * nmax);
                H5Dread(dsid, h5t_native<Integer>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
                count[i] = dset_data[0];
                free(dset_data);
            } // if
            H5Gclose(lvlid);
        } //  strcmp
    }     //  i
    H5Gclose(gid);
    H5Fclose(file);
    for (int i = count.size() - 1; i >= 0; --i) {
        if (count[i] == 0) {
            count.erase(count.begin() + i);
        }
    }
}

/// \brief          Get cells per level.
/// \brief          Get cells per level.
/// \tparam         Parameter Parameter type.
/// \tparam         Integer Count type
/// \tparam         Conic Cone type
/// \param[in]      parameters Parameter structure.
/// \param[in]      fileName File name.
/// \param[in,out]  count Number of cells per level
/// \param[in]      cubeNumber Name number of cubes per level.
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
/// \param[in]      conic Cone.
template <class Parameter, typename Integer, class Conic>
void TReadHDF5::cellsAndCubesPerLevels(const Parameter &parameters, const std::string &fileName, std::vector<Integer> &count, std::vector<std::vector<std::string>> &cubeNumber, const double &thetay, const double &thetaz, const Conic &conic) {

    float cubesize(0);
    TReadHDF5::getAttribute(fileName, "metadata/conecreator_grav_parameters/output_parameters", "cube_size", cubesize);
    std::vector<float> dimtab(3, 0);
    TReadHDF5::getAttribute(fileName, "metadata", "dimensions_array", dimtab);
    std::vector<int> nctab(3, 0);
    TReadHDF5::getAttribute(fileName, "metadata", "ncube_array", nctab);

    const double tanAngle = std::tan(conic.angle());
    const float normCone(conic.length());
    std::vector<float> baseCone(3, 0);
    std::vector<float> vertexCone(3, 0);
    baseCone[0] = conic.base(0);
    baseCone[1] = conic.base(1);
    baseCone[2] = conic.base(2);
    vertexCone[0] = conic.vertex(0);
    vertexCone[1] = conic.vertex(1);
    vertexCone[2] = conic.vertex(2);

    const double halfcubediag(0.5 * std::sqrt(3) * cubesize), microsphere(1. / parameters.microcoeff), coarsediag(std::sqrt(3) / (8 * parameters.microcoeff));
    const unsigned int MAX_NAME = 1024;
    char memb_name[MAX_NAME];
    H5G_info_t grpinfo, lvlinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, "/data/", H5P_DEFAULT);

    H5Gget_info(gid, &grpinfo);
    count.resize(grpinfo.nlinks);
    cubeNumber.resize(grpinfo.nlinks);
    // Loop over groups in '/data'
    for (unsigned int i = 0; i < grpinfo.nlinks; i++) {
        H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, memb_name, (size_t)MAX_NAME, H5P_DEFAULT);
        std::string levelName(memb_name);
        // Enter in groups which start with 'level'
        if (levelName.substr(0, 5) == "level") {
            hid_t lvlid = H5Gopen(gid, memb_name, H5P_DEFAULT);
            H5Gget_info(lvlid, &lvlinfo);
            if (lvlinfo.nlinks != 0) {
                cubeNumber.push_back(std::vector<std::string>());
                count.push_back(0);
                // Loop over levels
                for (unsigned int i1 = 0; i1 < lvlinfo.nlinks; i1++) {
                    H5Lget_name_by_idx(lvlid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name, (size_t)MAX_NAME, H5P_DEFAULT);
                    std::string str(memb_name);
                    // Select groups which start with 'cube'
                    if (str.substr(0, 4) == "cube") {
                        std::array<double, 3> point111;
                        if (parameters.isfullsky == 1) {
                            TReadHDF5::pos_from_index_fullsky(str, nctab, dimtab, cubesize, point111);
                        } else {
                            TReadHDF5::pos_from_index_narrow(str, nctab, dimtab, cubesize, thetay, thetaz, point111);
                        }
                        double length(0), distance(0), dist111(0);
                        for (unsigned int idim = 0; idim < 3; ++idim) {
                            length += (baseCone[idim] - vertexCone[idim]) * (point111[idim] - vertexCone[idim]);
                            dist111 += pow((point111[idim] - vertexCone[idim]), 2);
                        }
                        length = length / normCone;
                        for (unsigned int idim = 0; idim < 3; ++idim) {
                            distance += pow(point111[idim] - (vertexCone[idim] + (baseCone[idim] - vertexCone[idim]) * (length / normCone)), 2);
                        }
                        // Compute if the cube intersects the cone
                        if ((std::sqrt(dist111) < (halfcubediag + microsphere)) || ((std::sqrt(distance) < (5 * coarsediag + halfcubediag + length * tanAngle)) && (length > 0))) {
                            cubeNumber[i].push_back(str);
                            hid_t cubeid = H5Gopen(lvlid, memb_name, H5P_DEFAULT);
                            hid_t dsid = H5Dopen(cubeid, "refined_bool", H5P_DEFAULT);
                            hsize_t alloc = H5Dget_storage_size(dsid);
                            count[i] += alloc / (sizeof(Integer));
                            H5Dclose(dsid);
                            H5Gclose(cubeid);
                        } // if selection
                    }     // if cube
                }         // i1
            }             // if level
            H5Gclose(lvlid);
        } // strcmp
    }     //  i
    for (int i = count.size() - 1; i >= 0; --i) {
        // If no cells in any group at some given level, then erase level
        if (count[i] == 0) {
            count.erase(count.begin() + i);
            cubeNumber.erase(cubeNumber.begin() + i);
        }
    }
    H5Gclose(gid);
    H5Fclose(file);
}
// Attributs niveau 0
/// \brief          Get attribute.
/// \details        Get attribute.
/// \tparam         Type Attribute type.
/// \param[in]      fileName File name.
/// \param[in]	    attributeName Attribute name.
/// \param[in,out]  attributeValue Value of attribute.
template <typename Type>
void TReadHDF5::getAttribute(const std::string &fileName, const std::string &attributeName, Type &attributeValue) {

    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t aid = H5Aopen(file, attributeName.c_str(), H5P_DEFAULT);

    H5Aread(aid, h5t_native<Type>(), &attributeValue);

    H5Aclose(aid);
    H5Fclose(file);
}

/// \brief          Get attribute.
/// \details        Get attribute.
/// \tparam         Type Attribute type.
/// \param[in]      fileName File name.
/// \param[in]	    attributeName Attribute name.
/// \param[in,out]  attributeValue Values of attribute.
template <typename Type>
void TReadHDF5::getAttribute(const std::string &fileName, const std::string &attributeName, std::vector<Type> &attributeValue) {

    Type tabValues[3];
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t aid = H5Aopen(file, attributeName.c_str(), H5P_DEFAULT);

    H5Aread(aid, h5t_native<Type>(), &tabValues);

    attributeValue[0] = tabValues[0];
    attributeValue[1] = tabValues[1];
    attributeValue[2] = tabValues[2];

    H5Aclose(aid);
    H5Fclose(file);
}

// Attributs niveau 1
/// \brief          Get attribute.
/// \details        Get attribute.
/// \tparam         Type Attribute type.
/// \param[in]      fileName File name.
/// \param[in]	    groupName Group name.
/// \param[in] 	    attributeName Attribute name
/// \param[in,out]  attributeValue Value of attribute.
template <typename Type>
void TReadHDF5::getAttribute(const std::string &fileName, const std::string &groupName, const std::string &attributeName, Type &attributeValue) {

    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, groupName.c_str(), H5P_DEFAULT);
    hid_t aid = H5Aopen(gid, attributeName.c_str(), H5P_DEFAULT);

    H5Aread(aid, h5t_native<Type>(), &attributeValue);

    H5Aclose(aid);
    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get attribute.
/// \details        Get attribute.
/// \tparam         Type Attribute type.
/// \param[in]      fileName File name.
/// \param[in]      groupName Group name.
/// \param[in]	    attributeName Attribute name.
/// \param[in,out]  attributeValueValues of attribute.
template <typename Type>
void TReadHDF5::getAttribute(const std::string &fileName, const std::string &groupName, const std::string &attributeName, std::vector<Type> &attributeValue) {

    Type tabValues[3];
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, groupName.c_str(), H5P_DEFAULT);
    hid_t aid = H5Aopen(gid, attributeName.c_str(), H5P_DEFAULT);

    H5Aread(aid, h5t_native<Type>(), &tabValues);

    attributeValue[0] = tabValues[0];
    attributeValue[1] = tabValues[1];
    attributeValue[2] = tabValues[2];

    H5Aclose(aid);
    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from grav cones.
/// \details        Get list of values of a dataset.
/// \tparam         Type Data type.
/// \param[in]      fileName File name.
/// \param[in]      levelMin int Minimum level.
/// \param[in]      levelMax int Maximum level.
/// \param[in]	    output_name Data name.
/// \param[in,out]  output Data vector.
template <typename Type>
void TReadHDF5::fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::string &output_name, std::vector<Type> &output) {

    // Need to pre-compute the level names
    const std::vector<std::string> levels = {"level00", "level01", "level02", "level03", "level04", "level05", "level06", "level07", "level08", "level09", "level10", "level11", "level12", "level13", "level14", "level15", "level16", "level17", "level18", "level19", "level20", "level21", "level22", "level23", "level24"};
    unsigned int const MAX_NAME = 1024;
    char memb_name[MAX_NAME];
    H5G_info_t grpinfo, lvlinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    // Data must be in file at '/data'
    hid_t gid = H5Gopen(file, "/data/", H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    // Loop over levels
    for (unsigned int i = 0; i < levelMax - levelMin + 1; i++) {
        hid_t lvlid = H5Gopen(gid, levels[levelMin + i].c_str(), H5P_DEFAULT);
        H5Gget_info(lvlid, &lvlinfo);
        // Loop over groups in 'level'
        for (unsigned int i1 = 0; i1 < lvlinfo.nlinks; i1++) {
            H5Lget_name_by_idx(lvlid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select groups which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(lvlid, memb_name, H5P_DEFAULT);
                TReadHDF5::get_data_from_dataset(cubeid, output_name, output);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(lvlid);
    } // i
    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from grav cones.
/// \details        Get list of values of a dataset.
/// \tparam         Type Data type.
/// \param[in]      fileName File name.
/// \param[in]      levelMin int Minimum level.
/// \param[in]      levelMax int Maximum level.
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
/// \param[in]	    output_names Variadic Data names.
/// \param[in,out]  outputs Variadic Data vectors.
template <typename Type, typename... String, typename... Vector>
void TReadHDF5::fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs) {

    // Need to pre-compute the level names
    const std::vector<std::string> levels = {"level00", "level01", "level02", "level03", "level04", "level05", "level06", "level07", "level08", "level09", "level10", "level11", "level12", "level13", "level14", "level15", "level16", "level17", "level18", "level19", "level20", "level21", "level22", "level23", "level24"};
    unsigned int const MAX_NAME = 1024;
    char memb_name[MAX_NAME];
    H5G_info_t grpinfo, lvlinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    // Data must be in file at '/data'
    hid_t gid = H5Gopen(file, "/data/", H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    // Loop over levels
    for (unsigned int i = 0; i < levelMax - levelMin + 1; i++) {
        hid_t lvlid = H5Gopen(gid, levels[levelMin + i].c_str(), H5P_DEFAULT);
        H5Gget_info(lvlid, &lvlinfo);
        // Loop over groups in 'level'
        for (unsigned int i1 = 0; i1 < lvlinfo.nlinks; i1++) {
            H5Lget_name_by_idx(lvlid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select groups which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(lvlid, memb_name, H5P_DEFAULT);
                TReadHDF5::get_data_from_dataset(cubeid, output_name, output, output_names..., outputs...);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(lvlid);
    } // i
    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from grave cones.
/// \details        Get list of values of a dataset
///		    of a cube that intersects the cone
///		    GIVEN the list of cubes as inputs
/// \tparam         Type Data type.
/// \param[in]      fileName File name.
/// \param[in]      levelMin Minimum level.
/// \param[in]      levelMax Maximum level.
/// \param[in,out]  cubeNumber Name of cubes that intersects the cone.
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
template <typename Type>
void TReadHDF5::fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::vector<std::string> &cubeNumber, const std::string &output_name, std::vector<Type> &output) {

    const std::vector<std::string> levels = {"level00", "level01", "level02", "level03", "level04", "level05", "level06", "level07", "level08", "level09", "level10", "level11", "level12", "level13", "level14", "level15", "level16", "level17", "level18", "level19", "level20", "level21", "level22", "level23", "level24"};
    H5G_info_t grpinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, "/data/", H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    // Loop over levels
    for (unsigned int i = 0; i < levelMax - levelMin + 1; i++) {
        hid_t lvlid = H5Gopen(gid, levels[levelMin + i].c_str(), H5P_DEFAULT);
        // Loop over cube names
        for (unsigned int i1 = 0; i1 < cubeNumber.size(); i1++) {
            hid_t cubeid = H5Gopen(lvlid, cubeNumber[i1].c_str(), H5P_DEFAULT);
            TReadHDF5::get_data_from_dataset(cubeid, output_name, output);
            H5Gclose(cubeid);
        } // i1
        H5Gclose(lvlid);
    } // i
    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from grav cones.
/// \details        Get lists of values of datasets
///		    of a cube that intersects the cone
///		    GIVEN the list of cubes as inputs
/// \tparam         Type1 Data type.
/// \tparam         Type2 Data type.
/// \param[in]      fileName File name.
/// \param[in]      levelMin Minimum level.
/// \param[in]      levelMax Maximum level.
/// \param[in,out]  cubeNumber Name of cubes that intersects the cone.
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
/// \param[in]	    output_names Variadic Data names.
/// \param[in,out]  outputs Variadic Data vectors.
template <typename Type, typename... String, typename... Vector>
void TReadHDF5::fillVectors_grav(const std::string &fileName, const unsigned int &levelMin, const unsigned int &levelMax, const std::vector<std::string> &cubeNumber, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs) {

    const std::vector<std::string> levels = {"level00", "level01", "level02", "level03", "level04", "level05", "level06", "level07", "level08", "level09", "level10", "level11", "level12", "level13", "level14", "level15", "level16", "level17", "level18", "level19", "level20", "level21", "level22", "level23", "level24"};
    H5G_info_t grpinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, "/data/", H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    // Loop over levels
    for (unsigned int i = 0; i < levelMax - levelMin + 1; i++) {
        hid_t lvlid = H5Gopen(gid, levels[levelMin + i].c_str(), H5P_DEFAULT);
        // Loop over cube names
        for (unsigned int i1 = 0; i1 < cubeNumber.size(); i1++) {
            hid_t cubeid = H5Gopen(lvlid, cubeNumber[i1].c_str(), H5P_DEFAULT);
            TReadHDF5::get_data_from_dataset(cubeid, output_name, output, output_names..., outputs...);
            H5Gclose(cubeid);
        } // i1
        H5Gclose(lvlid);
    } // i
    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from part cones.
/// \details        Get lists of values of datasets
/// \tparam         Type dataset type.
/// \param[in]      fileName File name.
/// \param[in]      fileSide side name('data' or 'metadata' in Raygal simulation HDF5 files).
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
template <typename Type>
void TReadHDF5::fillVectors_part(const std::string &fileName, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output) {

    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);

    TReadHDF5::get_data_from_dataset(gid, output_name, output);

    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from part cones.
/// \details        Get lists of values of datasets
/// \tparam         Type1 Dataset type.
/// \tparam         Type2 Dataset type.
/// \param[in]      fileName File name.
/// \param[in]      fileSide side name('data' or 'metadata' in Raygal simulation HDF5 files).
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
/// \param[in]	    output_names Variadic Data names.
/// \param[in,out]  outputs Variadic Data vectors.
template <typename Type, typename... String, typename... Vector>
void TReadHDF5::fillVectors_part(const std::string &fileName, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs) {

    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
    // First dataset
    TReadHDF5::get_data_from_dataset(gid, output_name, output, output_names..., outputs...);

    H5Gclose(gid);
    H5Fclose(file);
}

/// \brief          Get data from part cones.
/// \details        Get lists of values of datasets
/// \tparam         Type1 Dataset type.
/// \param[in]      fraction Probability under which we take the data
/// \param[in]      fileName File name.
/// \param[in]      fileSide side name( data or metadata ).
/// \param[in]      output_name1 Data name.
/// \param[in,out]  output1 Data vector.
template <typename Type1>
void TReadHDF5::fillVectors_part(const double &fraction, const std::string &fileName, const std::string &fileSide, const std::string &output_name1, std::vector<Type1> &output1) {

    // Need for random selection
    if (fraction < 1) {
        unsigned int const MAX_NAME = 1024;
        char memb_name[MAX_NAME];
        std::uniform_real_distribution<double> distribution(0, 1);
        H5G_info_t grpinfo;
        hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
        H5Gget_info(gid, &grpinfo);
        // Loop over groups
        for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
            H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select groups which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                // First dataset (position or velocity)
                hid_t dsid = H5Dopen(cubeid, output_name1.c_str(), H5P_DEFAULT);
                hsize_t alloc = H5Dget_storage_size(dsid);
                unsigned int nmax = alloc / (sizeof(Type1));
                Type1 *dset_data1;
                dset_data1 = (Type1 *)malloc(sizeof(Type1) * nmax);
                H5Dread(dsid, h5t_native<Type1>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data1);
                std::vector<unsigned int> randomization(nmax / 3);
                std::vector<unsigned int> randomization_tmp;
                std::iota(std::begin(randomization), std::end(randomization), 0);
                // Randomize vector indexes and only keep a fraction
#ifdef GCCBELOW7
                std::experimental::sample(randomization.begin(), randomization.end(), std::back_inserter(randomization_tmp), static_cast<unsigned int>(fraction * nmax / 3), std::mt19937{std::random_device{}()});
#else
                std::sample(randomization.begin(), randomization.end(), std::back_inserter(randomization_tmp), static_cast<unsigned int>(fraction * nmax / 3), std::mt19937{std::random_device{}()});
#endif
                std::for_each(randomization_tmp.begin(), randomization_tmp.end(), [=, &output1](int i) { output1.insert(output1.end(), &dset_data1[3 * i], &dset_data1[3 * (i + 1)]); });
                free(dset_data1);
                H5Dclose(dsid);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(gid);
        H5Fclose(file);
        // Take full sample
    } else {
        unsigned int const MAX_NAME = 1024;
        char memb_name[MAX_NAME];
        H5G_info_t grpinfo;
        hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
        H5Gget_info(gid, &grpinfo);
        // Loop over groups
        for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
            H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select groups which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                TReadHDF5::get_data_from_dataset(cubeid, output_name1, output1);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(gid);
        H5Fclose(file);
    }
}

/// \brief          Get data from part cones.
/// \details        Get lists of values of datasets
/// \tparam         Type1 Dataset type.
/// \tparam         Type2 Dataset type.
/// \param[in]      fraction probability under which we take the data
/// \param[in]      fileName File name.
/// \param[in]      fileSide side name( data or metadata ).
/// \param[in]      output_name1 Data name.
/// \param[in,out]  output1 Data vector.
/// \param[in]	    output_name2 Variadic Data name.
/// \param[in,out]  output2 Variadic Data vector.
template <typename Type1, typename Type2>
void TReadHDF5::fillVectors_part(const double &fraction, const std::string &fileName, const std::string &fileSide, const std::string &output_name1, std::vector<Type1> &output1, const std::string &output_name2, std::vector<Type2> &output2) {

    // Need for random selection
    if (fraction < 1) {
        unsigned int const MAX_NAME = 1024;
        char memb_name[MAX_NAME];
        std::uniform_real_distribution<double> distribution(0, 1);
        H5G_info_t grpinfo;
        hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
        H5Gget_info(gid, &grpinfo);
        // Loop over groups
        for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
            H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select groups which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                // First dataset (position)
                hid_t dsid = H5Dopen(cubeid, output_name1.c_str(), H5P_DEFAULT);
                hsize_t alloc = H5Dget_storage_size(dsid);
                unsigned int nmax = alloc / (sizeof(Type1));
                std::vector<unsigned int> randomization(nmax / 3);
                std::vector<unsigned int> randomization_tmp;
                std::iota(std::begin(randomization), std::end(randomization), 0);
                // Randomize vector indexes and only keep a fraction
#ifdef GCCBELOW7
                std::experimental::sample(randomization.begin(), randomization.end(), std::back_inserter(randomization_tmp), static_cast<unsigned int>(fraction * nmax / 3), std::mt19937{std::random_device{}()});
#else
                std::sample(randomization.begin(), randomization.end(), std::back_inserter(randomization_tmp), static_cast<unsigned int>(fraction * nmax / 3), std::mt19937{std::random_device{}()});
#endif
                Type1 *dset_data1;
                dset_data1 = (Type1 *)malloc(sizeof(Type1) * nmax);
                H5Dread(dsid, h5t_native<Type1>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data1);
                std::for_each(randomization_tmp.begin(), randomization_tmp.end(), [=, &output1](int i) { output1.insert(output1.end(), &dset_data1[3 * i], &dset_data1[3 * (i + 1)]); });
                free(dset_data1);
                H5Dclose(dsid);
                // Second dataset (velocity)
                dsid = H5Dopen(cubeid, output_name2.c_str(), H5P_DEFAULT);
                alloc = H5Dget_storage_size(dsid);
                nmax = alloc / (sizeof(Type2));
                Type2 *dset_data2;
                dset_data2 = (Type2 *)malloc(sizeof(Type2) * nmax);
                H5Dread(dsid, h5t_native<Type2>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data2);
                std::for_each(randomization_tmp.begin(), randomization_tmp.end(), [=, &output2](int i) { output2.insert(output2.end(), &dset_data2[3 * i], &dset_data2[3 * (i + 1)]); });
                free(dset_data2);
                H5Dclose(dsid);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(gid);
        H5Fclose(file);
        // Take full sample
    } else {
        unsigned int const MAX_NAME = 1024;
        char memb_name[MAX_NAME];
        H5G_info_t grpinfo;
        hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
        H5Gget_info(gid, &grpinfo);
        // Loop over groups
        for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
            H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select groups which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                // First dataset
                TReadHDF5::get_data_from_dataset(cubeid, output_name1, output1);
                // Second dataset
                TReadHDF5::get_data_from_dataset(cubeid, output_name2, output2);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(gid);
        H5Fclose(file);
    }
}

/// \brief          Get data from part cones.
/// \details        Get lists of values of datasets
/// \tparam         Type1 Dataset type.
/// \tparam         Type2 Dataset type.
/// \tparam         Type3 Dataset type.
/// \param[in]      fraction probability under which we take the data
/// \param[in]      fileName File name.
/// \param[in]      fileSide side name( data or metadata ).
/// \param[in]      output_name1 Data name.
/// \param[in,out]  output1 Data vector.
/// \param[in]      output_name2 Data name.
/// \param[in,out]  output2 Data vector.
/// \param[in]      output_name3 Data name.
/// \param[in,out]  output3 Data vector.
template <typename Type1, typename Type2, typename Type3>
void TReadHDF5::fillVectors_part(const double &fraction, const std::string &fileName, const std::string &fileSide, const std::string &output_name1, std::vector<Type1> &output1, const std::string &output_name2, std::vector<Type2> &output2, const std::string &output_name3, std::vector<Type3> &output3) {

    // Need for random selection
    if (fraction < 1) {
        unsigned int const MAX_NAME = 1024;
        char memb_name[MAX_NAME];
        std::uniform_real_distribution<double> distribution(0, 1);
        H5G_info_t grpinfo;
        hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
        H5Gget_info(gid, &grpinfo);
        // Loop over groups
        for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
            H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select groups which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                // First dataset  (position)
                hid_t dsid = H5Dopen(cubeid, output_name1.c_str(), H5P_DEFAULT);
                hsize_t alloc = H5Dget_storage_size(dsid);
                unsigned int nmax1 = alloc / (sizeof(Type1));
                std::vector<unsigned int> randomization(nmax1 / 3);
                std::vector<unsigned int> randomization_tmp;
                std::iota(std::begin(randomization), std::end(randomization), 0);
                // Randomize vector indexes and only keep a fraction
#ifdef GCCBELOW7
                std::experimental::sample(randomization.begin(), randomization.end(), std::back_inserter(randomization_tmp), static_cast<unsigned int>(fraction * nmax1 / 3), std::mt19937{std::random_device{}()});
#else
                std::sample(randomization.begin(), randomization.end(), std::back_inserter(randomization_tmp), static_cast<unsigned int>(fraction * nmax1 / 3), std::mt19937{std::random_device{}()});
#endif
                Type1 *dset_data1;
                dset_data1 = (Type1 *)malloc(sizeof(Type1) * nmax1);
                H5Dread(dsid, h5t_native<Type1>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data1);
                std::for_each(randomization_tmp.begin(), randomization_tmp.end(), [=, &output1](int i) { output1.insert(output1.end(), &dset_data1[3 * i], &dset_data1[3 * (i + 1)]); });
                free(dset_data1);
                H5Dclose(dsid);
                // Second dataset (velocity)
                dsid = H5Dopen(cubeid, output_name2.c_str(), H5P_DEFAULT);
                alloc = H5Dget_storage_size(dsid);
                unsigned int nmax2 = alloc / (sizeof(Type2));
                Type2 *dset_data2;
                dset_data2 = (Type2 *)malloc(sizeof(Type2) * nmax2);
                H5Dread(dsid, h5t_native<Type2>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data2);
                std::for_each(randomization_tmp.begin(), randomization_tmp.end(), [=, &output2](int i) { output2.insert(output2.end(), &dset_data2[3 * i], &dset_data2[3 * (i + 1)]); });
                free(dset_data2);
                H5Dclose(dsid);
                // Third dataset (id)
                dsid = H5Dopen(cubeid, output_name3.c_str(), H5P_DEFAULT);
                alloc = H5Dget_storage_size(dsid);
                unsigned int nmax3 = alloc / (sizeof(Type3));
                Type3 *dset_data3;
                dset_data3 = (Type3 *)malloc(sizeof(Type3) * nmax3);
                H5Dread(dsid, h5t_native<Type3>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data3);
                std::for_each(randomization_tmp.begin(), randomization_tmp.end(), [=, &output3](int i) { output3.push_back(dset_data3[i]); });
                free(dset_data3);
                H5Dclose(dsid);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(gid);
        H5Fclose(file);
        // Take full sample
    } else {
        unsigned int const MAX_NAME = 1024;
        char memb_name[MAX_NAME];
        H5G_info_t grpinfo;
        hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
        H5Gget_info(gid, &grpinfo);
        // Loop over groups
        for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
            H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                               (size_t)MAX_NAME, H5P_DEFAULT);
            std::string str(memb_name);
            // Select group which start with 'cube'
            if (str.substr(0, 4) == "cube") {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                // First dataset
                TReadHDF5::get_data_from_dataset(cubeid, output_name1, output1);
                // Second dataset
                TReadHDF5::get_data_from_dataset(cubeid, output_name2, output2);
                // Third dataset
                TReadHDF5::get_data_from_dataset(cubeid, output_name3, output3);
                H5Gclose(cubeid);
            } // if cube
        }     // i1
        H5Gclose(gid);
        H5Fclose(file);
    }
}

/// \brief          Get data from part cones.
/// \details        Get values of datasets
/// \tparam         Parameter Parameter type.
/// \tparam         Type1 Dataset type.
/// \tparam         Type2 Dataset type.
/// \tparam         Conic Cone type.
/// \param[in]      parameters Parameter structure.
/// \param[in]      fileName File name.
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
/// \param[in]      conic Cone.
/// \param[in]      fileSide side name( data or metadata ).
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
template <class Parameter, typename Type, class Conic>
void TReadHDF5::fillVectors_part(const Parameter &parameters, const std::string &fileName, const double &thetay, const double &thetaz, const Conic &conic, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output) {

    float cubesize(0);
    TReadHDF5::getAttribute(fileName, "metadata/conecreator_part_parameters/output_parameters", "cube_size", cubesize);
    std::vector<float> dimtab(3, 0);
    TReadHDF5::getAttribute(fileName, "metadata", "dimensions_array", dimtab);
    std::vector<int> nctab(3, 0);
    TReadHDF5::getAttribute(fileName, "metadata", "ncube_array", nctab);

    const double tanAngle = std::tan(conic.angle());
    const float normCone(conic.length());
    std::vector<float> baseCone(3, 0);
    std::vector<float> vertexCone(3, 0);
    baseCone[0] = conic.base(0);
    baseCone[1] = conic.base(1);
    baseCone[2] = conic.base(2);
    vertexCone[0] = conic.vertex(0);
    vertexCone[1] = conic.vertex(1);
    vertexCone[2] = conic.vertex(2);

    double halfcubediag(0.5 * std::sqrt(3) * cubesize), microsphere(1. / parameters.microcoeff), coarsediag(std::sqrt(3) / (8 * parameters.microcoeff));
    unsigned int const MAX_NAME = 1024;
    char memb_name[MAX_NAME];
    H5G_info_t grpinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    // Loop over groups
    for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
        H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                           (size_t)MAX_NAME, H5P_DEFAULT);
        std::string str(memb_name);
        // Select groups which start with 'cube'
        if (str.substr(0, 4) == "cube") {
            std::array<double, 3> point111;
            // Different index calculation if fullsky or narrow cone
            if (parameters.isfullsky == 1) {
                TReadHDF5::pos_from_index_fullsky(str, nctab, dimtab, cubesize, point111);
            } else {
                TReadHDF5::pos_from_index_narrow(str, nctab, dimtab, cubesize, thetay, thetaz, point111);
            }
            double length(0), distance(0), dist111(0);
            for (unsigned int idim = 0; idim < 3; ++idim) {
                length += (baseCone[idim] - vertexCone[idim]) * (point111[idim] - vertexCone[idim]);
                dist111 += pow((point111[idim] - vertexCone[idim]), 2);
            }
            length = length / normCone;
            for (unsigned int idim = 0; idim < 3; ++idim) {
                distance += pow(point111[idim] - (vertexCone[idim] + (baseCone[idim] - vertexCone[idim]) * (length / normCone)), 2);
            }
            // Compute if cube intersects the cone
            if ((std::sqrt(dist111) < halfcubediag + microsphere) || ((std::sqrt(distance) < 2 * coarsediag + halfcubediag + length * tanAngle) && (length > 0))) {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                TReadHDF5::get_data_from_dataset(cubeid, output_name, output);
                H5Gclose(cubeid);
            } // if inside cone
        }     // if cube
    }         // i1
    H5Gclose(gid);
    H5Fclose(file);
}
/// \brief          Get data from part cones.
/// \details        Get values of datasets
/// \tparam         Parameter Parameter type.
/// \tparam         Type1 Dataset type.
/// \tparam         Type2 Dataset type.
/// \tparam         Conic Cone type.
/// \param[in]      parameters Parameter structure.
/// \param[in]      fileName File name.
/// \param[in]      thetay Semi-angle for solid angle in direction y
/// \param[in]      thetaz Semi-angle for solid angle in direction z
/// \param[in]      conic Cone.
/// \param[in]      fileSide side name( data or metadata ).
/// \param[in]      output_name Data name.
/// \param[in,out]  output Data vector.
/// \param[in]	    output_names Variadic Data names.
/// \param[in,out]  outputs Variadic Data vectors.
template <class Parameter, typename Type, typename... String, typename... Vector, class Conic>
void TReadHDF5::fillVectors_part(const Parameter &parameters, const std::string &fileName, const double &thetay, const double &thetaz, const Conic &conic, const std::string &fileSide, const std::string &output_name, std::vector<Type> &output, const String &...output_names, Vector &...outputs) {

    float cubesize(0);
    TReadHDF5::getAttribute(fileName, "metadata/conecreator_part_parameters/output_parameters", "cube_size", cubesize);
    std::vector<float> dimtab(3, 0);
    TReadHDF5::getAttribute(fileName, "metadata", "dimensions_array", dimtab);
    std::vector<int> nctab(3, 0);
    TReadHDF5::getAttribute(fileName, "metadata", "ncube_array", nctab);

    const double tanAngle = std::tan(conic.angle());
    const float normCone(conic.length());
    std::vector<float> baseCone(3, 0);
    std::vector<float> vertexCone(3, 0);
    baseCone[0] = conic.base(0);
    baseCone[1] = conic.base(1);
    baseCone[2] = conic.base(2);
    vertexCone[0] = conic.vertex(0);
    vertexCone[1] = conic.vertex(1);
    vertexCone[2] = conic.vertex(2);

    double halfcubediag(0.5 * std::sqrt(3) * cubesize), microsphere(1. / parameters.microcoeff), coarsediag(std::sqrt(3) / (8 * parameters.microcoeff));
    unsigned int const MAX_NAME = 1024;
    char memb_name[MAX_NAME];
    H5G_info_t grpinfo;
    hid_t file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen(file, fileSide.c_str(), H5P_DEFAULT);
    H5Gget_info(gid, &grpinfo);
    // Loop over groups
    for (unsigned int i1 = 0; i1 < grpinfo.nlinks; i1++) {
        H5Lget_name_by_idx(gid, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i1, memb_name,
                           (size_t)MAX_NAME, H5P_DEFAULT);
        std::string str(memb_name);
        // Select groups which start with 'cube'
        if (str.substr(0, 4) == "cube") {
            std::array<double, 3> point111;
            // Different index calculation if fullsky or narrow cone
            if (parameters.isfullsky == 1) {
                TReadHDF5::pos_from_index_fullsky(str, nctab, dimtab, cubesize, point111);
            } else {
                TReadHDF5::pos_from_index_narrow(str, nctab, dimtab, cubesize, thetay, thetaz, point111);
            }
            double length(0), distance(0), dist111(0);
            for (unsigned int idim = 0; idim < 3; ++idim) {
                length += (baseCone[idim] - vertexCone[idim]) * (point111[idim] - vertexCone[idim]);
                dist111 += pow((point111[idim] - vertexCone[idim]), 2);
            }
            length = length / normCone;
            for (unsigned int idim = 0; idim < 3; ++idim) {
                distance += pow(point111[idim] - (vertexCone[idim] + (baseCone[idim] - vertexCone[idim]) * (length / normCone)), 2);
            }
            // Compute if cube intersects the cone
            if ((std::sqrt(dist111) < halfcubediag + microsphere) || ((std::sqrt(distance) < 2 * coarsediag + halfcubediag + length * tanAngle) && (length > 0))) {
                hid_t cubeid = H5Gopen(gid, memb_name, H5P_DEFAULT);
                TReadHDF5::get_data_from_dataset(cubeid, output_name, output, output_names..., outputs...);
                H5Gclose(cubeid);
            } // if inside cone
        }     // if cube
    }         // i1
    H5Gclose(gid);
    H5Fclose(file);
}

#endif // TReadHDF5_H_INCLUDED
