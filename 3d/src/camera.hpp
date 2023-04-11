#pragma once

#include "mve/math/math.hpp"
#include "mve/window.hpp"

class Camera {
public:
    Camera();

    void update(const mve::Window& window);

    void fixed_update(const mve::Window& window);

    inline mve::Matrix4 view_matrix(float interpolation_weight) const
    {
        mve::Matrix4 transform = m_head_transform * m_body_transform;
        mve::Matrix3 basis = transform.basis();
        mve::Matrix4 interpolated_transform = mve::Matrix4::from_basis_translation(
            basis, m_prev_pos.linear_interpolate(transform.translation(), interpolation_weight));
        mve::Matrix4 view = interpolated_transform.inverse().transpose();
        return view;
    }

private:
    mve::Matrix4 m_body_transform;
    mve::Matrix4 m_head_transform;
    mve::Vector3 m_prev_pos;
    float m_friction;
    float m_acceleration;
    float m_max_speed;
    mve::Vector3 m_velocity;
};