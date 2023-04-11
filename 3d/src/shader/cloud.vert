#version 460

layout (set = 0, binding = 0) uniform GlobalUniform {
    mat4 view;
    mat4 proj;
} global_ubo;

layout (location = 0) in vec3 in_pos;

void main() {
    gl_Position = global_ubo.proj * global_ubo.view * vec4(in_pos, 1.0);
}