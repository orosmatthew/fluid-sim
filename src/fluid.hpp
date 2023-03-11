#pragma once

#include <algorithm>
#include <cassert>
#include <functional>
#include <vector>

class Fluid {
public:
    inline Fluid(int size, float diffusion, float viscosity, int iter)
    {
        assert(size > 0);
        assert(diffusion >= 0.0f);
        assert(viscosity >= 0.0f);
        assert(iter > 0);

        m_size = size;
        m_diff = diffusion;
        m_viscosity = viscosity;

        std::vector<float> empty(size * size, 0.0f);
        m_s = empty;
        m_density = empty;

        m_vel_x = empty;
        m_vel_y = empty;

        m_vel_x_next = empty;
        m_vel_y_next = empty;

        m_pressure = empty;
        m_divergence = empty;

        m_lin_solve_iterations = iter;
    }

    inline void step(float time_step)
    {
        for (size_t i = 0; i < m_size * m_size; i++) {
            m_density[i] = std::clamp(m_density[i], 0.0f, 10000.0f);
        }

        diffuse(BoundaryType::neumann, m_vel_x, m_vel_x_next, m_viscosity, time_step, m_size, m_lin_solve_iterations);
        diffuse(BoundaryType::neumann, m_vel_y, m_vel_y_next, m_viscosity, time_step, m_size, m_lin_solve_iterations);

        calc_pressure(m_vel_x_next, m_vel_y_next, m_pressure, m_divergence, m_size, m_lin_solve_iterations);
        correct_velocity(m_vel_x_next, m_vel_y_next, m_pressure, m_size);

        advect(BoundaryType::neumann, m_vel_x_next, m_vel_x, m_vel_x_next, m_vel_y_next, time_step, m_size);
        advect(BoundaryType::neumann, m_vel_y_next, m_vel_y, m_vel_x_next, m_vel_y_next, time_step, m_size);

        calc_pressure(m_vel_x, m_vel_y, m_pressure, m_divergence, m_size, m_lin_solve_iterations);
        correct_velocity(m_vel_x, m_vel_y, m_pressure, m_size);

        diffuse(BoundaryType::fixed, m_density, m_s, m_diff, time_step, m_size, m_lin_solve_iterations);
        advect(BoundaryType::fixed, m_s, m_density, m_vel_x, m_vel_y, time_step, m_size);

        for (size_t i = 0; i < m_size * m_size; i++) {
            m_density[i] = std::clamp(m_density[i], 0.0f, 10000.0f);
        }
    }

    inline void add_velocity(int x, int y, float amount_x, float amount_y)
    {
        size_t i = index(x, y, m_size);
        m_vel_x[i] += amount_x;
        m_vel_y[i] += amount_y;
    }

    inline void add_density(int x, int y, float amount)
    {
        m_density[index(x, y, m_size)] += amount;
    }

    [[nodiscard]] inline float density_at(int x, int y) const
    {
        return m_density[index(x, y, m_size)];
    }

    [[nodiscard]] inline float pressure_at(int x, int y) const
    {
        return m_pressure[index(x, y, m_size)];
    }

    [[nodiscard]] inline int size() const
    {
        return m_size;
    }

private:
    struct Vector2 {
        float x;
        float y;
    };

    struct Vector2i {
        int x;
        int y;
    };

    int m_size;
    float m_diff;
    float m_viscosity;

    std::vector<float> m_s;
    std::vector<float> m_density;

    std::vector<float> m_vel_x;
    std::vector<float> m_vel_y;

    std::vector<float> m_vel_x_next;
    std::vector<float> m_vel_y_next;

    std::vector<float> m_pressure;
    std::vector<float> m_divergence;

    int m_lin_solve_iterations;

    inline static float lerp(float a, float b, float t)
    {
        return a * (1 - t) + b * t;
    };

    inline static void for_2d(const Vector2i& min, const Vector2i& max, std::function<void(Vector2i)> func)
    {
        for (int i = min.x; i < max.x; i++) {
            for (int j = min.y; j < max.y; j++) {
                std::invoke(func, Vector2i { i, j });
            }
        }
    }

    inline static size_t index(int x, int y, int size)
    {
        return static_cast<size_t>(x) + static_cast<size_t>(y) * size;
    }

    enum class BoundaryType {
        none, // No boundary condition, used for free surface boundaries
        fixed, // Fixed boundary condition, aka Dirichlet condition
        neumann // Neumann boundary condition, aka zero-gradient condition
    };

    inline static void set_bnd(BoundaryType boundary_type, std::vector<float>& x, int size)
    {
        // Boundary conditions for left and right
        for (int i = 1; i < size - 1; i++) {
            if (boundary_type == BoundaryType::neumann) { // Set boundary values to negative of adjacent value
                x[index(i, 0, size)] = -x[index(i, 1, size)];
                x[index(i, size - 1, size)] = -x[index(i, size - 2, size)];
            }
            else if (boundary_type == BoundaryType::fixed) { // Set boundary values to same value as the adjacent value
                x[index(i, 0, size)] = x[index(i, 1, size)];
                x[index(i, size - 1, size)] = x[index(i, size - 2, size)];
            }
        }
        // Boundary conditions for top and bottom
        for (int j = 1; j < size - 1; j++) {
            if (boundary_type == BoundaryType::neumann) {
                x[index(0, j, size)] = -x[index(1, j, size)];
                x[index(size - 1, j, size)] = -x[index(size - 2, j, size)];
            }
            else if (boundary_type == BoundaryType::fixed) {
                x[index(0, j, size)] = x[index(1, j, size)];
                x[index(size - 1, j, size)] = x[index(size - 2, j, size)];
            }
        }
        // Boundary condition for corners
        if (boundary_type == BoundaryType::fixed) {
            x[index(0, 0, size)] = 0.5f * (x[index(1, 0, size)] + x[index(0, 1, size)]);
            x[index(0, size - 1, size)] = 0.5f * (x[index(1, size - 1, size)] + x[index(0, size - 2, size)]);
            x[index(size - 1, 0, size)] = 0.5f * (x[index(size - 2, 0, size)] + x[index(size - 1, 1, size)]);
            x[index(size - 1, size - 1, size)]
                = 0.5f * (x[index(size - 2, size - 1, size)] + x[index(size - 1, size - 2, size)]);
        }
        // No boundary condition for corners is needed for Neumann or free surfaces
    }

    // Solve for scalar field for Poisson equation.
    // The Poisson equation is a partial differential equation that relates the second derivative of a scalar
    // function to a source term.
    // L^2*u=f
    // L^2 = Laplacian operator; u = scalar field (dest); f = scalar field (src)
    // The Laplacian operator describes the rate at which a scalar field changes over space (u -> f)
    // This function essentially computes the previous scalar field given the current one
    inline static void lin_solve(
        BoundaryType boundary_type,
        std::vector<float>& dest,
        const std::vector<float>& src,
        float a,
        float c,
        int size,
        int iter)
    {
        const float c_inv = 1.0f / c;

        for (int t = 0; t < iter; t++) {
            for_2d({ 1, 1 }, { size - 1, size - 1 }, [&](const Vector2i& pos) {
                // clang-format off
                float neighbor_sum =
                    dest[index(pos.x + 1, pos.y, size)] +
                    dest[index(pos.x - 1, pos.y, size)] +
                    dest[index(pos.x, pos.y + 1, size)] +
                    dest[index(pos.x, pos.y - 1, size)];
                // clang-format on

                // Contribution of the laplacian operator to the new value of the current point
                float laplacian_contribution = a * neighbor_sum;

                float new_value = (src[index(pos.x, pos.y, size)] + laplacian_contribution) * c_inv;

                dest[index(pos.x, pos.y, size)] = new_value;
            });
            set_bnd(boundary_type, dest, size);
        }
    }

    inline static void calc_pressure(
        const std::vector<float>& vel_x,
        const std::vector<float>& vel_y,
        std::vector<float>& pressure,
        std::vector<float>& divergence,
        int size,
        int iter)
    {
        // Calculate the divergence of the velocity field.
        // Divergence in the velocity field is a scalar value that measures how much the fluid is flowing outward or
        // inward at a given point.
        // (+) away from point
        // (-) toward point
        // (0) not accumulating nor depleting
        // This is used to compute the pressure field which helps to correct for numerical errors by conserving
        // mass, momentum and energy.
        for_2d({ 1, 1 }, { size - 1, size - 1 }, [&](const Vector2i pos) {
            const float delta_x = vel_x[index(pos.x + 1, pos.y, size)] - vel_x[index(pos.x - 1, pos.y, size)];
            const float delta_y = vel_y[index(pos.x, pos.y + 1, size)] - vel_y[index(pos.x, pos.y - 1, size)];
            const size_t current = index(pos.x, pos.y, size);
            divergence[current] = -0.5f * (delta_x + delta_y) / static_cast<float>(size);
            pressure[current] = 0;
        });

        set_bnd(BoundaryType::none, divergence, size);
        set_bnd(BoundaryType::none, pressure, size);

        // Scalar constant is used to scale the result of the Poisson equation before updating pressure
        const float scalar_constant = 1.5f;

        // The discretization constant is needed to discretize continuous equations.
        // Represents spacing or step size between points
        const float discretization_constant = 6.0f;

        // Solve for pressure by solving the Poisson equation using the divergence of the velocity field as the source.
        lin_solve(BoundaryType::none, pressure, divergence, scalar_constant, discretization_constant, size, iter);
    }

    inline static void correct_velocity(
        std::vector<float>& vel_x, std::vector<float>& vel_y, const std::vector<float>& pressure, int size)
    {
        // Subtract the pressure gradient from the velocity field to ensure incompressibility
        // This ensures that the velocity field remains divergence-free
        for_2d({ 1, 1 }, { size - 1, size - 1 }, [&](const Vector2i& pos) {
            const size_t current = index(pos.x, pos.y, size);
            vel_x[current] -= 0.5f * (pressure[index(pos.x + 1, pos.y, size)] - pressure[index(pos.x - 1, pos.y, size)])
                * static_cast<float>(size);
            vel_y[current] -= 0.5f * (pressure[index(pos.x, pos.y + 1, size)] - pressure[index(pos.x, pos.y - 1, size)])
                * static_cast<float>(size);
        });

        set_bnd(BoundaryType::neumann, vel_x, size);
        set_bnd(BoundaryType::neumann, vel_y, size);
    }

    inline static void advect(
        BoundaryType boundary_type,
        const std::vector<float>& from,
        std::vector<float>& to,
        const std::vector<float>& vel_x,
        const std::vector<float>& vel_y,
        float time_step,
        int size)
    {
        Vector2 dt { time_step * static_cast<float>(size) - 2, time_step * static_cast<float>(size) - 2 };

        for_2d({ 1, 1 }, { size - 1, size - 1 }, [&](const Vector2i& pos) {
            // Index of current pos
            const size_t current = index(pos.x, pos.y, size);
            // displacement = dt * vel
            const Vector2 displacement { dt.x * vel_x[current], dt.y * vel_y[current] };
            // new_pos = pos - displacement
            Vector2 new_pos { static_cast<float>(pos.x) - displacement.x, static_cast<float>(pos.y) - displacement.y };
            // Clamp new position to size
            new_pos.x = std::clamp(new_pos.x, 0.5f, static_cast<float>(size) - 2 + 0.5f);
            new_pos.y = std::clamp(new_pos.y, 0.5f, static_cast<float>(size) - 2 + 0.5f);
            // new_pos_i = int(floor(new_pos))
            const Vector2i new_pos_i { static_cast<int>(floorf(new_pos.x)), static_cast<int>(floorf(new_pos.y)) };
            // offset = new_pos - new_pos_i
            const Vector2 offset { new_pos.x - static_cast<float>(new_pos_i.x),
                                   new_pos.y - static_cast<float>(new_pos_i.y) };

            // clang-format off
            // Neighboring indices of points in the direction of the displacement/velocity
            const size_t neighbors[2][2] = {
                {index(new_pos_i.x,     new_pos_i.y, size), index(new_pos_i.x,     new_pos_i.y + 1, size)},
                {index(new_pos_i.x + 1, new_pos_i.y, size), index(new_pos_i.x + 1, new_pos_i.y + 1, size)}
            };

            // Perform bilinear interpolation between neighbors
            to[current] =
                (1.0f - offset.x) * lerp(from[neighbors[0][0]], from[neighbors[0][1]], offset.y) +
                offset.x          * lerp(from[neighbors[1][0]], from[neighbors[1][1]], offset.y);
            // clang-format on
        });
        set_bnd(boundary_type, to, size);
    }

    inline static void diffuse(
        BoundaryType boundary_type,
        const std::vector<float>& from,
        std::vector<float>& to,
        float diffusion_constant,
        float time_step,
        int size,
        int iter)
    {
        // Scaling factor for linear solving.
        // Higher value results in faster rate of diffusion.
        float a = time_step * diffusion_constant * (static_cast<float>(size) - 2) * (static_cast<float>(size) - 2);
        lin_solve(boundary_type, to, from, a, 1 + 6 * a, size, iter);
    }
};