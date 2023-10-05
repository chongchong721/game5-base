#include "Mode.hpp"

#include "Scene.hpp"
#include "WalkMesh.hpp"

#include <glm/glm.hpp>

#include <vector>
#include <deque>

#include <chrono>

struct PlayMode : Mode {
	PlayMode();
	virtual ~PlayMode();

	//functions called by main loop:
	virtual bool handle_event(SDL_Event const &, glm::uvec2 const &window_size) override;
	virtual void update(float elapsed) override;
	virtual void draw(glm::uvec2 const &drawable_size) override;

	void return_to_walkmesh();
	void init_player();

	//----- game state -----

	//input tracking:
	struct Button {
		uint8_t downs = 0;
		uint8_t pressed = 0;
	} left, right, speed_down, speed_up;

	//local copy of the game scene (so code can change it during gameplay):
	Scene scene;

	//player info:
	struct Player {
		WalkPoint at;
		//transform is at player's feet and will be yawed by mouse left/right motion:
		Scene::Transform *transform = nullptr;
		//camera is at player's head and will be pitched by mouse up/down motion:
		Scene::Camera *camera = nullptr;
		glm::vec3 speed{-0.01f,0.0f,0.0f};
		bool on_walkmesh = true;
		float max_speed = 2.0f;
		glm::vec3 car_direction = {-1.0f,0.0f,0.0f};
		Scene::Transform *original_transform = nullptr;
		uint8_t reset_time = 3;
	} player;

	struct{
		glm::vec3 min;
		glm::vec3 max;
	}bbox;

	glm::vec3 gravity{0.0,0.0,-0.5};
	bool failed = false;
	bool success = false;
	float exit_timer = 5.0f;
	std::chrono::time_point<std::chrono::system_clock> * start = nullptr;
	std::chrono::time_point<std::chrono::system_clock> end;
	std::chrono::duration<double> elapsed_time;

};
