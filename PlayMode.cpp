#include "PlayMode.hpp"

#include "LitColorTextureProgram.hpp"
#include "frameProgram.hpp"

#include "DrawLines.hpp"
#include "Mesh.hpp"
#include "Load.hpp"
#include "gl_errors.hpp"
#include "data_path.hpp"

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>

#include <random>

GLuint meshes_for_lit_color_texture_program = 0;
Load< MeshBuffer > meshes(LoadTagDefault, []() -> MeshBuffer const * {
	MeshBuffer const *ret = new MeshBuffer(data_path("mesh.pnct"));
	meshes_for_lit_color_texture_program = ret->make_vao_for_program(lit_color_texture_program->program);
	return ret;
});

Load< Scene > mesh_scene(LoadTagDefault, []() -> Scene const * {
	return new Scene(data_path("mesh.scene"), [&](Scene &scene, Scene::Transform *transform, std::string const &mesh_name){
		Mesh const &mesh = meshes->lookup(mesh_name);

		scene.drawables.emplace_back(new Scene::Drawable(transform));
		std::shared_ptr<Scene::Drawable> drawable = scene.drawables.back();

		drawable->pipeline = lit_color_texture_program_pipeline;

		drawable->pipeline.vao = meshes_for_lit_color_texture_program;
		drawable->pipeline.type = mesh.type;
		drawable->pipeline.start = mesh.start;
		drawable->pipeline.count = mesh.count;

	});
});

WalkMesh const *walkmesh = nullptr;
Load< WalkMeshes > walkmeshes(LoadTagDefault, []() -> WalkMeshes const * {
	WalkMeshes *ret = new WalkMeshes(data_path("mesh.w"));
	walkmesh = &ret->lookup("Mesh");
	return ret;
});

PlayMode::PlayMode() : scene(*mesh_scene) {
	//create a player transform:
	//scene.transforms.emplace_back();
	for(auto &t : scene.transforms){
		if (t.name == "car"){
			player.transform = &t;
			player.original_transform = new Scene::Transform();
			*player.original_transform = t;
		}
	}


	for(auto &d : scene.drawables){
		if(d->transform->name == "car"){
			d->draw_frame = true;
			// auto pipeline = d->pipeline;
			// d->pipeline = frame_program_pipeline;
			// d->pipeline.vao = pipeline.vao;
			// d->pipeline.count = pipeline.count;
			// d->pipeline.type = pipeline.type;
			// d->pipeline.start = pipeline.start;
			// d->pipeline.count = pipeline.count;
			// for (auto i  = 0U ; i < pipeline.TextureCount; i++){
			// 	d->pipeline.textures[i] = pipeline.textures[i];
			// }
		}
	}

	auto m = meshes->lookup("end");
	bbox.min = m.min;
	bbox.max = m.max;

	//create a player camera attached to a child of the player transform:
	scene.transforms.emplace_back();
	scene.cameras.emplace_back(&scene.transforms.back());
	player.camera = &scene.cameras.back();
	player.camera->fovy = glm::radians(60.0f);
	player.camera->near = 0.01f;
	player.camera->transform->parent = player.transform;

	//debug
	//player.transform->position = glm::vec3(0.6,8.7,1.9);

	//player's eyes are 0.5 units above the ground and behind the car:
	player.camera->transform->position = glm::vec3(0.4f, 0.0f, 0.15f);

	//rotate camera facing direction (-z) to player facing direction (-x):
	player.camera->transform->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
	player.camera->transform->rotation *= glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));


	//start player walking at nearest walk point:
	player.at = walkmesh->nearest_walk_point(player.transform->position);

}

PlayMode::~PlayMode() {
}

bool PlayMode::handle_event(SDL_Event const &evt, glm::uvec2 const &window_size) {

	if (evt.type == SDL_KEYDOWN) {
		if (evt.key.keysym.sym == SDLK_ESCAPE) {
			SDL_SetRelativeMouseMode(SDL_FALSE);
			return true;
		} else if (evt.key.keysym.sym == SDLK_a) {
			left.downs += 1;
			left.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.downs += 1;
			right.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			speed_up.downs += 1;
			speed_up.pressed = true;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			speed_down.downs += 1;
			speed_down.pressed = true;
			return true;
		}
	} else if (evt.type == SDL_KEYUP) {
		if (evt.key.keysym.sym == SDLK_a) {
			left.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_d) {
			right.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_w) {
			speed_up.pressed = false;
			return true;
		} else if (evt.key.keysym.sym == SDLK_s) {
			speed_down.pressed = false;
			return true;
		}
	} else if (evt.type == SDL_MOUSEBUTTONDOWN) {
		if (SDL_GetRelativeMouseMode() == SDL_FALSE) {
			SDL_SetRelativeMouseMode(SDL_TRUE);
			return true;
		}
	} else if (evt.type == SDL_MOUSEMOTION) {
		// if (SDL_GetRelativeMouseMode() == SDL_TRUE) {
		// 	glm::vec2 motion = glm::vec2(
		// 		evt.motion.xrel / float(window_size.y),
		// 		-evt.motion.yrel / float(window_size.y)
		// 	);
		// 	glm::vec3 upDir = walkmesh->to_world_smooth_normal(player.at);
		// 	player.transform->rotation = glm::angleAxis(-motion.x * player.camera->fovy, upDir) * player.transform->rotation;

		// 	float pitch = glm::pitch(player.camera->transform->rotation);
		// 	pitch += motion.y * player.camera->fovy;
		// 	//camera looks down -z (basically at the player's feet) when pitch is at zero.
		// 	pitch = std::min(pitch, 0.95f * 3.1415926f);
		// 	pitch = std::max(pitch, 0.05f * 3.1415926f);
		// 	player.camera->transform->rotation = glm::angleAxis(pitch, glm::vec3(1.0f, 0.0f, 0.0f));

		// 	return true;
		// }
		return true;
	}

	return false;
}

void PlayMode::update(float elapsed) {

	if(success || failed){
		end = std::chrono::system_clock::now();
		elapsed_time = end - *start;
		if(elapsed_time.count() > exit_timer){
			Mode::set_current(nullptr);
		}
	} 
	

	// Assume the car has a constant accelration

	if(player.on_walkmesh)
	//player walking:
	{
		glm::vec3 curreint_pos = player.transform->position;
		glm::vec3 move = curreint_pos;

		// If reach the destination?
		if (curreint_pos[0] > bbox.min[0] && curreint_pos[1]> bbox.min[1] && curreint_pos[0] < bbox.max[0] && curreint_pos[1] < bbox.max[1] ){
			if(!success){
				success = true;
				if(start == nullptr)
				{	
					start = new std::chrono::time_point<std::chrono::system_clock>;
					*start = std::chrono::system_clock::now();
				}
			}
		}





		//rotate the car? compute the acceleration
		glm::vec3 acc{0.0f};	
		glm::vec3 dir = player.car_direction;
		auto normal = walkmesh->to_world_triangle_normal(player.at);
		if(left.pressed && !right.pressed){
			auto rotation = glm::angleAxis(glm::radians(0.2f),normal);
			dir = rotation * glm::vec4(dir,0.0);
			player.car_direction = dir;
			player.transform->rotation *= rotation;
		}
		if(!left.pressed && right.pressed){
			auto rotation = glm::angleAxis(glm::radians(-0.2f),normal);
			dir = rotation * glm::vec4(dir,0.0);
			player.car_direction = dir;
			player.transform->rotation *= rotation;
		}

		//combine inputs into a move:
		//constexpr float PlayerSpeed = 3.0f;





		// gravity along car direction
		glm::vec3 g = gravity - (glm::dot(gravity,normal) * normal);

		if (speed_down.pressed && !speed_up.pressed){
			acc = -dir * glm::vec3(0.2) + g;
			move = curreint_pos + player.speed * elapsed + glm::vec3(0.5f) * acc * elapsed * elapsed;
			player.speed += acc * elapsed;
		} else if  (!speed_down.pressed && speed_up.pressed){
			acc = dir * glm::vec3(0.2) + g;
			move = curreint_pos + player.speed * elapsed + glm::vec3(0.5f) * acc * elapsed * elapsed;
			player.speed += acc * elapsed;

		} else{
			acc = g;
			move = curreint_pos + player.speed * elapsed + glm::vec3(0.5f) * acc * elapsed * elapsed;
			player.speed += acc * elapsed;
		}
		if(glm::length(player.speed) > player.max_speed){
			player.speed *= glm::vec3{player.max_speed / glm::length(player.speed)};
		}


		//get move in world coordinate system:
		//glm::vec3 remain = player.transform->make_local_to_world() * glm::vec4(move.x, move.y, move.z, 0.0f);
		glm::vec3 remain = move - curreint_pos;

		//using a for() instead of a while() here so that if walkpoint gets stuck in
		// some awkward case, code will not infinite loop:
		for (uint32_t iter = 0; iter < 10; ++iter) {
			if (remain == glm::vec3(0.0f)) break;
			WalkPoint end;
			float time;
			walkmesh->walk_in_triangle(player.at, remain, &end, &time);
			player.at = end;
			if (time == 1.0f) {
				//finished within triangle:
				remain = glm::vec3(0.0f);
				break;
			}
			//some step remains:
			remain *= (1.0f - time);
			//try to step over edge:
			glm::quat rotation;
			if (walkmesh->cross_edge(player.at, &end, &rotation)) {
				//stepped to a new triangle:
				player.at = end;
				//rotate step to follow surface:
				remain = rotation * remain;
				//rotate speed
				player.speed = rotation * glm::vec4(player.speed,0.0f);
				// rotate direction
				player.car_direction = rotation * player.car_direction;
			} else {
				//TODO: just let it free-fall.
				//ran into a wall, bounce / slide along it:
				glm::vec3 const &a = walkmesh->vertices[player.at.indices.x];
				glm::vec3 const &b = walkmesh->vertices[player.at.indices.y];
				glm::vec3 const &c = walkmesh->vertices[player.at.indices.z];
				glm::vec3 along = glm::normalize(b-a);
				glm::vec3 normal = glm::normalize(glm::cross(b-a, c-a));
				glm::vec3 in = glm::cross(normal, along);

				//check how much 'remain' is pointing out of the triangle:
				float d = glm::dot(remain, in);
				if (d < 0.0f) {
					//Go out?
					player.on_walkmesh = false;
				} else {
					//if it's just pointing along the edge, bend slightly away from wall:
					remain += 0.01f * d * in;
				}
			}

			player.speed = player.car_direction * glm::dot(player.car_direction,player.speed);
		}

		if (remain != glm::vec3(0.0f)) {
			std::cout << "NOTE: code used full iteration budget for walking." << std::endl;
		}

		//update player's position to respect walking:
		player.transform->position = walkmesh->to_world_point(player.at);

		{ //update player's rotation to respect local (smooth) up-vector:
			
			glm::quat adjust = glm::rotation(
				player.transform->rotation * glm::vec3(0.0f, 0.0f, 1.0f), //current up vector
				walkmesh->to_world_smooth_normal(player.at) //smoothed up vector at walk location
			);
			player.transform->rotation = glm::normalize(adjust * player.transform->rotation);
		}

		/*
		glm::mat4x3 frame = camera->transform->make_local_to_parent();
		glm::vec3 right = frame[0];
		//glm::vec3 up = frame[1];
		glm::vec3 forward = -frame[2];

		camera->transform->position += move.x * right + move.y * forward;
		*/
	}else{
		glm::vec3 curreint_pos = player.transform->position;
		if (curreint_pos[2] <= -5.0f){
			init_player();
		}else{
			glm::vec3 acc = gravity;
			glm::vec3 move = curreint_pos + player.speed * elapsed + glm::vec3(0.5f) * acc * elapsed * elapsed;
			move = move - curreint_pos;
			player.speed += acc * elapsed;
			// Rotate the player?

			// How to rotate to make the car perpendicular to Z-axis?
			glm::vec3 dir = glm::vec3{0.0,0.0,1.0} * glm::dot(player.car_direction,glm::vec3(0.0,0.0,1.0));
			if(dir.z != 0.0){
				glm::quat adjust = glm::rotation(
					player.transform->rotation * glm::vec3(0.0f, 0.0f, 1.0f), //current up vector
					glm::vec3(0.0,0.0,1.0) //smoothed up vector at walk location
				);
				player.transform->rotation = glm::normalize(adjust * player.transform->rotation);
				player.car_direction = player.transform->rotation * glm::vec4{player.car_direction,0.0};
				player.car_direction = glm::normalize(glm::vec3(player.car_direction.x,player.car_direction.y,0.0));
			}
			player.transform->position = player.transform->position + move;

			//Find the nearest walkmesh.
			return_to_walkmesh();
		}
		


	}

	//reset button press counters:
	left.downs = 0;
	right.downs = 0;
	speed_up.downs = 0;
	speed_down.downs = 0;
}

void PlayMode::draw(glm::uvec2 const &drawable_size) {
	//update camera aspect ratio for drawable:
	player.camera->aspect = float(drawable_size.x) / float(drawable_size.y);

	//set up light type and position for lit_color_texture_program:
	// TODO: consider using the Light(s) in the scene to do this
	glUseProgram(lit_color_texture_program->program);
	glUniform1i(lit_color_texture_program->LIGHT_TYPE_int, 1);
	glUniform3fv(lit_color_texture_program->LIGHT_DIRECTION_vec3, 1, glm::value_ptr(glm::vec3(0.0f, 0.0f,-1.0f)));
	glUniform3fv(lit_color_texture_program->LIGHT_ENERGY_vec3, 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 0.95f)));
	glUseProgram(0);

	glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
	glClearDepth(1.0f); //1.0 is actually the default value to clear the depth buffer to, but FYI you can change it.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS); //this is the default depth comparison function, but FYI you can change it.


	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	scene.draw(*player.camera, false);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	scene.draw(*player.camera, true);
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);

	/* In case you are wondering if your walkmesh is lining up with your scene, try:
	{
		glDisable(GL_DEPTH_TEST);
		DrawLines lines(player.camera->make_projection() * glm::mat4(player.camera->transform->make_world_to_local()));
		for (auto const &tri : walkmesh->triangles) {
			lines.draw(walkmesh->vertices[tri.x], walkmesh->vertices[tri.y], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
			lines.draw(walkmesh->vertices[tri.y], walkmesh->vertices[tri.z], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
			lines.draw(walkmesh->vertices[tri.z], walkmesh->vertices[tri.x], glm::u8vec4(0x88, 0x00, 0xff, 0xff));
		}
	}
	*/

	{ //use DrawLines to overlay some text:
		glDisable(GL_DEPTH_TEST);
		float aspect = float(drawable_size.x) / float(drawable_size.y);
		DrawLines lines(glm::mat4(
			1.0f / aspect, 0.0f, 0.0f, 0.0f,
			0.0f, 1.0f, 0.0f, 0.0f,
			0.0f, 0.0f, 1.0f, 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f
		));

		std::string text;

		if(!success && !failed){
			text = "WS moves; AD turns";
		}else if (success){
			char out[100];
			auto length = sprintf(out,"Congrats! You Win! Game will end in %.1f seconds",exit_timer - elapsed_time.count());
			text = std::string(out,length);
		}else if (failed){
			char out[100];
			auto length = sprintf(out,"You failed! Game will end in %.1f seconds",exit_timer - elapsed_time.count());
			text = std::string(out,length);
		}

		constexpr float H = 0.09f;
		lines.draw_text(text,
			glm::vec3(-aspect + 0.1f * H, -1.0 + 0.1f * H, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0x00, 0x00, 0x00, 0x00));
		float ofs = 2.0f / drawable_size.y;
		lines.draw_text(text,
			glm::vec3(-aspect + 0.1f * H + ofs, -1.0 + + 0.1f * H + ofs, 0.0),
			glm::vec3(H, 0.0f, 0.0f), glm::vec3(0.0f, H, 0.0f),
			glm::u8vec4(0xff, 0xff, 0xff, 0x00));


	}
	GL_ERRORS();
}


void PlayMode::return_to_walkmesh(){
	glm::vec3 curreint_pos = player.transform->position;
	auto wp = walkmesh->nearest_walk_point(curreint_pos);
	
	auto in_01 = [](float x)->bool{
		if(x > 0.0 && x < 1.0)
			return true;
		else
			return false;
	};

	if(in_01(wp.weights[0]) && in_01(wp.weights[1]) && in_01(wp.weights[2])){
		auto pt = walkmesh->to_world_point(wp);
		//printf("z1,z0:%.3f,%.3f\n",curreint_pos[2],pt[2]);
		if ((curreint_pos[2] - pt[2]) >= 0 && (curreint_pos[2]-pt[2]) < 0.1){
			player.on_walkmesh = true;
			player.transform->position = player.transform->position + (pt - curreint_pos);
			player.at = wp;
			// Let the car keep moving
			player.speed[2] = 0.0f;
			//printf("return to walkmesh\n");
		}
	}
}


void PlayMode::init_player(){
	if(!success && player.reset_time == 0){
		failed = true;
		if(start == nullptr){
			
			start = new std::chrono::time_point<std::chrono::system_clock>;
			*start = std::chrono::system_clock::now();
				
		}
	}

	*player.transform = *player.original_transform;
	player.camera->fovy = glm::radians(60.0f);
	player.camera->near = 0.01f;
	player.camera->transform->parent = player.transform;

	//player's eyes are 0.5 units above the ground and behind the car:
	player.camera->transform->position = glm::vec3(0.4f, 0.0f, 0.15f);

	//rotate camera facing direction (-z) to player facing direction (-x):
	player.camera->transform->rotation = glm::angleAxis(glm::radians(90.0f), glm::vec3(0.0f, 0.0f, 1.0f));
	player.camera->transform->rotation *= glm::angleAxis(glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));


	//start player walking at nearest walk point:
	player.at = walkmesh->nearest_walk_point(player.transform->position);
	player.speed = glm::vec3{-0.01f,0.0f,0.0f};
	player.car_direction = glm::vec3{-1.0f,0.0f,0.0f};
	player.on_walkmesh = true;
	player.reset_time -= 1;
}