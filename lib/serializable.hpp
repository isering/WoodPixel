/*
WoodPixel - Supplementary code for Computational Parquetry:
            Fabricated Style Transfer with Wood Pixels
            ACM Transactions on Graphics 39(2), 2020

Copyright (C) 2020 Julian Iseringhausen <opensource@iseringhausen.graphics>
Copyright (C) 2020 Matthias Hullin, University of Bonn <hullin@cs.uni-bonn.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef TRLIB_SERIALIZABLE_HPP_
#define TRLIB_SERIALIZABLE_HPP_

#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <opencv2/opencv.hpp>

class Serializable
{
public:
  virtual void load(const boost::filesystem::path& base_path, const boost::property_tree::ptree& tree) = 0;
  virtual boost::property_tree::ptree save(const boost::filesystem::path& base_path, const boost::filesystem::path& path) const = 0;
  
  static boost::filesystem::path get_unique_path(const boost::filesystem::path& base_path, const boost::filesystem::path& path, const std::string& extension)
  {
    if (!boost::filesystem::exists(base_path / path))
    {
      boost::filesystem::create_directories(base_path / path);
    }

    int i = 0;
    boost::filesystem::path path_unique;
    boost::filesystem::path path_relative;
    do
    {
      std::stringstream ss;
      ss << std::setw(8) << std::setfill('0') << i++ << extension;
      path_unique = base_path / path / ss.str();
      path_relative = path / ss.str();
    } while (boost::filesystem::exists(path_unique));

    return path_relative;
  }

  /*
  * std::vector of Serializable
  */
  template <typename T, typename std::enable_if<std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, const std::vector<T>& vec, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_vec;
    for (typename std::vector<T>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter)
    {
      tree_vec.push_back(std::make_pair("", iter->save(base_path, path)));
    }
    tree.add_child(key, tree_vec);
  }

  template <typename Iter, typename std::enable_if<std::is_base_of<Serializable, typename Iter::value_type>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, Iter begin, Iter end, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_vec;
    for (Iter iter = begin; iter != end; ++iter)
    {
      tree_vec.push_back(std::make_pair("", iter->save(base_path, path)));

    }
    tree.add_child(key, tree_vec);
  }

  template <typename T, typename std::enable_if<std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, std::vector<T>& vec, const boost::filesystem::path& base_path)
  {
    vec.clear();
    for (const auto& tree_val : tree.get_child(key))
    {
      T val;
      deserialize(tree_val.second, "", val, base_path);
      vec.push_back(val);
    }
  }

  /*
  * std::vector of cv::Point
  */
  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, const std::vector<cv::Point_<T>>& vec, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_vec;
    for (typename std::vector<cv::Point_<T>>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter)
    {
      boost::property_tree::ptree tree_elem;
      serialize(tree_elem, "x", iter->x, base_path, path);
      serialize(tree_elem, "y", iter->y, base_path, path);
      tree_vec.push_back(std::make_pair("", tree_elem));
    }
    tree.add_child(key, tree_vec);
  }

  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, std::vector<cv::Point_<T>>& vec, const boost::filesystem::path& base_path)
  {
    vec.clear();
    for (const auto& val : tree.get_child(key))
    {
      vec.emplace_back(val.second.get<T>("x"), val.second.get<T>("y"));
    }
  }

  /*
  * std::vector of Generic
  */
  template <typename Iter, typename std::enable_if<!std::is_base_of<Serializable, typename Iter::value_type>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, Iter begin, Iter end, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_vec;
    for (Iter iter = begin; iter != end; ++iter)
    {
      boost::property_tree::ptree tree_elem;
      tree_elem.put("", *iter);
      tree_vec.push_back(std::make_pair("", tree_elem));
    }
    tree.add_child(key, tree_vec);
  }

  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, std::vector<T>& vec, const boost::filesystem::path& base_path)
  {
    vec.clear();
    for (const auto& val : tree.get_child(key))
    {
      vec.push_back(val.second.get_value<T>());
    }
  }

  /*
   * Serializable
   */
  template <typename T, typename std::enable_if<std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, const T& child, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    tree.add_child(key, child.save(base_path, path));
  }

  template <typename T, typename std::enable_if<std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, T& child, const boost::filesystem::path& base_path)
  {
    child.load(base_path, tree.get_child(key));
  }

  /*
   * cv::Point
   */
  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, const cv::Point_<T>& p, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_point;
    serialize(tree_point, "x", p.x, base_path, path);
    serialize(tree_point, "y", p.y, base_path, path);
    tree.add_child(key, tree_point);
  }

  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, cv::Point_<T>& p, const boost::filesystem::path& base_path)
  {
    const boost::property_tree::ptree tree_point = tree.get_child(key);
    deserialize(tree_point, "x", p.x, base_path);
    deserialize(tree_point, "y", p.y, base_path);
  }

  /*
   * cv::Rect
   */
  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, const cv::Rect_<T>& rect, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_rect;
    serialize(tree_rect, "x", rect.x, base_path, path);
    serialize(tree_rect, "y", rect.y, base_path, path);
    serialize(tree_rect, "width", rect.width, base_path, path);
    serialize(tree_rect, "height", rect.height, base_path, path);
    tree.add_child(key, tree_rect);
  }

  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, cv::Rect_<T>& rect, const boost::filesystem::path& base_path)
  {
    const boost::property_tree::ptree tree_rect = tree.get_child(key);
    deserialize(tree_rect, "x", rect.x, base_path);
    deserialize(tree_rect, "y", rect.y, base_path);
    deserialize(tree_rect, "width", rect.width, base_path);
    deserialize(tree_rect, "height", rect.height, base_path);
  }

  /*
   * cv::Size
   */
  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, const cv::Size_<T>& size, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_size;
    serialize(tree_size, "width", size.width, base_path, path);
    serialize(tree_size, "height", size.height, base_path, path);
    tree.add_child(key, tree_size);
  }

  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, cv::Size_<T>& size, const boost::filesystem::path& base_path)
  {
    const boost::property_tree::ptree tree_size = tree.get_child(key);
    deserialize(tree_size, "width", size.width, base_path);
    deserialize(tree_size, "height", size.height, base_path);
  }

  /*
   * Generic
   */
  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void serialize(boost::property_tree::ptree& tree, const std::string& key, const T& val, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    tree.add(key, val);
  }

  template <typename T, typename std::enable_if<!std::is_base_of<Serializable, T>::value>::type* = nullptr>
  static void deserialize(const boost::property_tree::ptree& tree, const std::string& key, T& val, const boost::filesystem::path& base_path)
  {
    val = tree.get<T>(key);
  }

  /*
   * Enums
   */
  template <typename T>
  static void serialize_enum(boost::property_tree::ptree& tree, const std::string& key, const T& val, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    tree.add(key, static_cast<int>(val));
  }

  template <typename T>
  static void deserialize_enum(const boost::property_tree::ptree& tree, const std::string& key, T& val, const boost::filesystem::path& base_path)
  {
    val = static_cast<T>(tree.get<int>(key));
  }

  /*
   * cv::Mat as image
   */
  static void serialize_image(boost::property_tree::ptree& tree, const std::string& key, const cv::Mat& mat, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    if (mat.empty())
    {
      serialize(tree, key, "<empty>", base_path, path);
    }
    else
    {
      const boost::filesystem::path path_tiff = get_unique_path(base_path, path, ".tiff");
      cv::imwrite((base_path / path_tiff).string(), mat);
      serialize(tree, key, path_tiff.string(), base_path, path);
    }
  }

  static void serialize_image_path(boost::property_tree::ptree& tree, const std::string& key, const boost::filesystem::path& path_image, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    cv::Mat image = cv::imread(path_image.string(), cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH);
    serialize_image(tree, key, image, base_path, path);
  }

  static void deserialize_image(const boost::property_tree::ptree& tree, const std::string& key, cv::Mat& mat, const boost::filesystem::path& base_path)
  {
    std::string filename;
    deserialize(tree, key, filename, base_path);
    if (filename == "<empty>")
    {
      mat = cv::Mat();
    }
    else
    {
      mat = cv::imread((base_path / filename).string(), cv::IMREAD_ANYDEPTH | cv::IMREAD_ANYCOLOR);
    }
  }

  /*
   * cv::Mat 
  */
  template <typename T>
  static void serialize_mat(boost::property_tree::ptree& tree, const std::string& key, const cv::Mat& mat, const boost::filesystem::path& base_path, const boost::filesystem::path& path)
  {
    boost::property_tree::ptree tree_mat;
    serialize(tree_mat, "rows", mat.rows, base_path, path);
    serialize(tree_mat, "cols", mat.cols, base_path, path);
    serialize(tree_mat, "channels", mat.channels(), base_path, path);
    serialize(tree_mat, "mat", mat.begin<T>(), mat.end<T>(), base_path, path);
    tree.add_child(key, tree_mat);
  }

  template <typename T>
  static void deserialize_mat(const boost::property_tree::ptree& tree, const std::string& key, cv::Mat& mat, const boost::filesystem::path& base_path)
  {
    std::vector<T> data;
    int channels, rows;
    
    const boost::property_tree::ptree tree_mat = tree.get_child(key);
    deserialize(tree_mat, "mat", data, base_path);
    deserialize(tree_mat, "channels", channels, base_path);
    deserialize(tree_mat, "rows", rows, base_path);
    
    mat = cv::Mat_<T>(data, true);
    mat = mat.reshape(channels, rows);
  }

  template <typename T>
  static std::vector<T> deserialize_vec(const boost::property_tree::ptree& tree, const std::string& key, const boost::filesystem::path& base_path)
  {
    std::vector<T> vec;
    for (const auto& tree_patch : tree.get_child(key))
    {
      vec.emplace_back();
      vec.back().load(base_path, tree_patch.second);
    }
    return vec;
  }
};

#endif /* TRLIB_SERIALIZABLE_HPP_ */