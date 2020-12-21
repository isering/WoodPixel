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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "texture.hpp"

#include <stack>
#include <stdexcept>
#include <regex>

#include <boost/filesystem.hpp>
#include <boost/math/constants/constants.hpp>

#include "affine_transformation.hpp"
#include "bezier_curve.hpp"
#include "bezier_ransac.hpp"
#include "hough_circle_detect.hpp"

namespace fs = boost::filesystem;

void Texture::load_texture(const boost::filesystem::path &filename, double scale)
{
  std::ifstream file(filename.string(), std::iostream::binary);
  if (!file.good())
  {
    throw std::runtime_error("Unable to read image: " + filename.string());
  }

  file.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
  file.seekg(0, std::ios::end);
  std::streampos length(file.tellg());
  std::vector<char> buffer(static_cast<std::size_t>(length));
  if (static_cast<std::size_t>(length) == 0)
  {
    throw std::runtime_error("Image file empty: " + filename.string());
  }
  file.seekg(0, std::ios::beg);
  file.read(buffer.data(), static_cast<std::size_t>(length));
  file.close();

  texture = cv::imdecode(buffer, cv::IMREAD_ANYDEPTH | cv::IMREAD_COLOR);

  if (!texture.data)
  {
    throw std::invalid_argument("Unable to load image: " + filename.string());
  }

  if (texture.depth() == CV_8U)
  {
    texture.convertTo(texture, CV_16UC3, 255.0);
  }
  else if (texture.depth() == CV_32F)
  {
    texture.convertTo(texture, CV_16UC3, 65535.0);
  }

  if (texture.depth() != CV_16U || texture.channels() != 3)
  {
    throw std::runtime_error("Texture loading failed");
  }
  cv::resize(texture, texture, cv::Size(), scale, scale);
}

void Texture::load_mask(const boost::filesystem::path &filename, double scale)
{
  mask_done = cv::imread(filename.string(), cv::IMREAD_GRAYSCALE);
  if (!mask_done.data)
  {
    throw(std::runtime_error("Unable to load mask: " + filename.string()));
  }
  for (auto iter = mask_done.begin<unsigned char>(); iter != mask_done.end<unsigned char>(); ++iter)
  {
    *iter = *iter > 0 ? 255 : 0;
  }
  cv::resize(mask_done, mask_done, cv::Size(), scale, scale, cv::INTER_NEAREST);
}

static cv::Rect2f compute_rotated_rect(const cv::Mat &texture, float angle_rad)
{
  const float angle_deg = 180.0f * angle_rad / boost::math::constants::pi<float>();
  return cv::RotatedRect(cv::Point2f(), texture.size(), angle_deg).boundingRect2f();
}

cv::Mat Texture::compute_transformation_matrix(const cv::Mat &texture, double angle_rad)
{
  const cv::Rect2f rotated_rect = compute_rotated_rect(texture, static_cast<float>(angle_rad));
  const double angle_deg = 180.0 * angle_rad / boost::math::constants::pi<double>();
  const cv::Point2f texture_center(0.5f * (texture.cols - 1), 0.5f * (texture.rows - 1));

  cv::Mat transformation_matrix = cv::getRotationMatrix2D(texture_center, angle_deg, 1.0);
  transformation_matrix.at<double>(0, 2) += 0.5 * rotated_rect.width - 0.5 * texture.cols;
  transformation_matrix.at<double>(1, 2) += 0.5 * rotated_rect.height - 0.5 * texture.rows;

  return transformation_matrix;
}

Texture Texture::rotate(double angle_rad) const
{
  cv::Rect2f rotated_rect = compute_rotated_rect(texture, static_cast<float>(angle_rad));

  Texture texture_rotated;
  texture_rotated.angle_rad = angle_rad;
  texture_rotated.dpi = dpi;
  texture_rotated.scale = scale;
  texture_rotated.transformation_matrix = compute_transformation_matrix(texture, angle_rad);
  cv::invertAffineTransform(texture_rotated.transformation_matrix, texture_rotated.transformation_matrix_inv);

  cv::warpAffine(texture, texture_rotated.texture, texture_rotated.transformation_matrix, rotated_rect.size());
  cv::warpAffine(mask_done, texture_rotated.mask_done, texture_rotated.transformation_matrix, rotated_rect.size(), cv::INTER_NEAREST);
  cv::warpAffine(mask_rotation, texture_rotated.mask_rotation, texture_rotated.transformation_matrix, rotated_rect.size(), cv::INTER_NEAREST);

  for (const cv::Point2d &p : marker.markers_pix)
  {
    const cv::Point2d p_transformed = AffineTransformation::transform(texture_rotated.transformation_matrix, p);
    texture_rotated.marker.markers_pix.push_back(p);
  }

  texture_rotated.filename = filename;

  return texture_rotated;
}

void Texture::mask_patch(std::vector<Texture> &rotated_textures, int index, cv::Point anchor, const std::vector<cv::Point> &patch)
{
  cv::Mat mask(rotated_textures[index].mask_done.size(), CV_8UC1, cv::Scalar(255));
  for (const cv::Point &p : patch)
  {
    mask.at<unsigned char>(anchor + p) = 0;
  }

  for (Texture &texture : rotated_textures)
  {
    cv::Mat transform_patch = AffineTransformation::concat(texture.transformation_matrix, rotated_textures[index].transformation_matrix_inv);
    cv::Mat mask_rotated;
    cv::warpAffine(mask, mask_rotated, transform_patch, texture.mask_done.size(), cv::INTER_NEAREST);
    mask_rotated = cv::max(mask_rotated, texture.mask_rotation_inv());
    texture.mask_done = cv::min(texture.mask_done, mask_rotated);
  }
}

Texture Texture::operator()(const cv::Range &row_range, const cv::Range &col_range)
{
  Texture texture_out;

  texture_out.texture = texture(row_range, col_range);
  texture_out.mask_done = mask_done(row_range, col_range);
  texture_out.mask_rotation = mask_rotation(row_range, col_range);
  texture_out.transformation_matrix = transformation_matrix;
  texture_out.transformation_matrix_inv = transformation_matrix_inv;
  texture_out.angle_rad = angle_rad;
  texture_out.scale = scale;
  texture_out.response = response(row_range, col_range);
  texture_out.filename = filename;

  return texture_out;
}

Texture Texture::operator()(const cv::Rect &rect)
{
  return (*this)(cv::Range(rect.y, rect.y + rect.height), cv::Range(rect.x, rect.x + rect.width));
}

/*
cv::Mat Texture::template_match_gpu(const Texture& kernel) const
{
  cv::Mat response_float, kernel_response_float;

  Timer<boost::milli> t;
  
  cv::cuda::Stream stream;

  const int rows_out = response.rows() - kernel.response.rows() + 1;
  const int cols_out = response.cols() - kernel.response.cols() + 1;

  std::vector<cv::cuda::GpuMat> texture_gpu(response.num_channels(), cv::cuda::GpuMat(response.size(), CV_32FC1));
  std::vector<cv::cuda::GpuMat> kernel_gpu(response.num_channels(), cv::cuda::GpuMat(response.size(), CV_32FC1));
  std::vector<cv::cuda::GpuMat> match_gpu(response.num_channels(), cv::cuda::GpuMat(response.size(), CV_32FC1));

  std::vector<cv::Ptr<cv::cuda::TemplateMatching>> matcher;

  boost::chrono::duration<float, boost::milli> duration_convert(0.0f);
  boost::chrono::duration<float, boost::milli> duration_upload(0.0f);
  boost::chrono::duration<float, boost::milli> duration_match(0.0f);

  for (int i = 0; i < response.num_channels(); ++i)
  {
    Timer<boost::milli> t_convert;
    response[i].convertTo(response_float, CV_32FC1, 1.0 / 65535.0);
    kernel.response[i].convertTo(kernel_response_float, CV_32FC1, 1.0 / 65535.0);
    duration_convert += t_convert.duration();

    Timer<boost::milli> t_upload;
    texture_gpu[i].upload(response_float);
    kernel_gpu[i].upload(kernel_response_float);
    duration_upload += t_upload.duration();

    Timer<boost::milli> t_match;
    matcher.push_back(cv::cuda::createTemplateMatching(CV_32F, CV_TM_SQDIFF));
    matcher[i]->match(texture_gpu[i], kernel_gpu[i], match_gpu[i]);
    duration_match += t_match.duration();
  }

  for (int i = 1; i < response.num_channels(); ++i)
  {
    cv::cuda::add(match_gpu[0], match_gpu[i], match_gpu[0], cv::noArray(), -1, stream);
  }

  cv::Mat match;
  match_gpu[0].download(match, stream);

  stream.waitForCompletion();

  t.stop();
  std::cout << "Duration template match GPU: " << t << " (" << duration_convert << " / " << duration_upload << " / " << duration_match << ")" << std::endl;

  return match;
}
*/

cv::Mat Texture::template_match(const Texture &kernel) const
{
  cv::Mat response_float, kernel_response_float;

  const int rows_out = response.rows() - kernel.response.rows() + 1;
  const int cols_out = response.cols() - kernel.response.cols() + 1;

  cv::Mat match(rows_out, cols_out, CV_32FC1);
  cv::Mat match_sum = cv::Mat::zeros(rows_out, cols_out, CV_32FC1);

  for (int i = 0; i < response.num_channels(); ++i)
  {
    response[i].convertTo(response_float, CV_32FC1, 1.0 / 65535.0);
    kernel.response[i].convertTo(kernel_response_float, CV_32FC1, 1.0 / 65535.0);
    cv::matchTemplate(response_float, kernel_response_float, match, cv::TM_SQDIFF);
    match_sum += match;
  }

  return match_sum;
}

cv::Mat Texture::template_match(const Texture &kernel, cv::Mat mask) const
{
  cv::Mat response_float, kernel_response_float;

  const int rows_out = response.rows() - kernel.response.rows() + 1;
  const int cols_out = response.cols() - kernel.response.cols() + 1;

  cv::Mat match(rows_out, cols_out, CV_32FC1);
  cv::Mat match_sum = cv::Mat::zeros(rows_out, cols_out, CV_32FC1);

  cv::Mat mask_float = cv::Mat::zeros(mask.size(), CV_32FC1);
  mask_float.setTo(1.0f, mask != 0);

  for (int i = 0; i < response.num_channels(); ++i)
  {
    response[i].convertTo(response_float, CV_32FC1, 1.0 / 65535.0);
    kernel.response[i].convertTo(kernel_response_float, CV_32FC1, 1.0 / 65535.0);
    cv::matchTemplate(response_float, kernel_response_float, match, cv::TM_SQDIFF, mask_float);
    match_sum += match;
  }

  return match_sum;
}

Texture Texture::clone() const
{
  Texture rhs;
  rhs.texture = texture.clone();
  rhs.mask_done = mask_done.clone();
  rhs.mask_rotation = mask_rotation.clone();
  rhs.transformation_matrix = transformation_matrix.clone();
  rhs.transformation_matrix_inv = transformation_matrix_inv.clone();
  rhs.angle_rad = angle_rad;
  rhs.dpi = dpi;
  rhs.scale = scale;
  rhs.response = response;
  rhs.filename = filename;
  return rhs;
}

void Texture::downsample_nn(int factor)
{
  double f = 1.0 / factor;
  cv::resize(texture, texture, cv::Size(), f, f, cv::INTER_AREA);
  cv::resize(mask_done, mask_done, cv::Size(), f, f, cv::INTER_NEAREST);
  cv::resize(mask_rotation, mask_rotation, cv::Size(), f, f, cv::INTER_NEAREST);
  response.downsample_nn(factor);
  scale /= factor;
}

std::vector<cv::Vec3f> Texture::find_markers(double marker_size_mm, int num_markers)
{
  std::vector<cv::Vec3f> markers;

  cv::Mat texture_8bit, texture_gray;
  texture.convertTo(texture_8bit, CV_8UC1, 1.0 / 255.0);
  cv::cvtColor(texture_8bit, texture_gray, cv::COLOR_BGR2GRAY);

  HoughCircleDetect circle_detector(0.5);
  circle_detector.compute_interactive(markers, texture_gray, texture_8bit, dpi, marker_size_mm, "Marker segmentation");

  return markers;
}

TextureRegion Texture::get_regions(cv::Rect region, cv::Mat edge_image)
{
  std::vector<cv::Point2d> edge_points;
  for (int y = region.y; y < region.y + region.height; ++y)
  {
    const unsigned char *ptr = reinterpret_cast<const unsigned char *>(edge_image.ptr(y));
    for (int x = region.x; x < region.x + region.width; ++x)
    {
      if (ptr[x])
      {
        edge_points.emplace_back(x, y);
      }
    }
  }

  if (edge_points.empty())
  {
    return TextureRegion();
  }

  BezierCurve curve = fit_bezier_cubic_ransac(edge_points, region, 6, 100);
  if (curve.is_empty_curve())
  {
    return TextureRegion();
  }

  cv::Mat mask = curve.compute_separating_curve_mask(region);

  /*
  mask.setTo(64, mask==1);
  mask.setTo(128, mask==2);
  cv::resize(mask, mask, cv::Size(), 8.0, 8.0, cv::INTER_NEAREST);
  cv::imshow("Mask", mask);
  cv::waitKey();
  std::exit(EXIT_SUCCESS);
  */

  const int num_mask_1 = cv::countNonZero(mask == 1);
  const int num_mask_2 = cv::countNonZero(mask == 2);

  TextureRegion region_out;
  region_out.mask_1 = (num_mask_1 <= num_mask_2 ? ((mask == 1) + (mask == 255)) : (mask == 1));
  region_out.mask_2 = (num_mask_1 > num_mask_2 ? ((mask == 2) + (mask == 255)) : (mask == 2));
  region_out.region = region;
  region_out.separating_curve = curve;

  return region_out;
}

cv::Vec3b Texture::interpolate_texture(const cv::Point2f &p) const
{
  const int x = static_cast<int>(p.x);
  const int y = static_cast<int>(p.y);
  const float dx = p.x - x;
  const float dy = p.y - y;

  const cv::Vec3f c11 = cv::Vec3f(texture.at<cv::Vec3w>(y, x)) / 65535.0f;
  const cv::Vec3f c12 = cv::Vec3f(texture.at<cv::Vec3w>(y, x + 1)) / 65535.0f;
  const cv::Vec3f c21 = cv::Vec3f(texture.at<cv::Vec3w>(y + 1, x)) / 65535.0f;
  const cv::Vec3f c22 = cv::Vec3f(texture.at<cv::Vec3w>(y + 1, x + 1)) / 65535.0f;
  const cv::Vec3f c1 = (1.0f - dx) * c11 + dx * c12;
  const cv::Vec3f c2 = (1.0f - dx) * c21 + dx * c22;
  const cv::Vec3f c = (1.0f - dy) * c1 + dy * c2;

  return cv::Vec3b(255.0f * c);
}

boost::property_tree::ptree Texture::save(const boost::filesystem::path &base_path, const boost::filesystem::path &path) const
{
  boost::property_tree::ptree tree;
  serialize_image(tree, "texture", texture, base_path, path);
  serialize_image(tree, "mask", mask_done, base_path, path);
  serialize_image(tree, "mask_rotation", mask_rotation, base_path, path);
  serialize_mat<double>(tree, "transformation_matrix", transformation_matrix, base_path, path);
  serialize_mat<double>(tree, "transformation_matrix_inv", transformation_matrix_inv, base_path, path);
  serialize(tree, "angle_rad", angle_rad, base_path, path);
  serialize(tree, "scale", scale, base_path, path);
  serialize(tree, "dpi", dpi, base_path, path);
  serialize(tree, "marker", marker, base_path, path);
  serialize(tree, "id", id, base_path, path);
  serialize_image_path(tree, "filename", filename, base_path, path);
  return tree;
}

void Texture::load(const boost::filesystem::path &base_path, const boost::property_tree::ptree &tree)
{
  deserialize_image(tree, "texture", texture, base_path);
  deserialize_image(tree, "mask", mask_done, base_path);
  deserialize_image(tree, "mask_rotation", mask_rotation, base_path);
  deserialize_mat<double>(tree, "transformation_matrix", transformation_matrix, base_path);
  deserialize_mat<double>(tree, "transformation_matrix_inv", transformation_matrix_inv, base_path);
  deserialize(tree, "angle_rad", angle_rad, base_path);
  deserialize(tree, "scale", scale, base_path);
  deserialize(tree, "dpi", dpi, base_path);
  deserialize(tree, "marker", marker, base_path);
  deserialize(tree, "id", id, base_path);
  deserialize(tree, "filename", filename, base_path);
}
