---
permalink: /instructors
layout: page
title: Instructors
---
  <section class="px-4">
  {% for instructor in site.instructors %}
	<figure class="image is-128x128" style="float:left; clear:left" >
		<img class="is-rounded" src="{{ site.baseurl }}/assets/instructors/{{ instructor.picture }}">
	</figure>
    <div class="px-4" style="overflow: auto">
    	<h2>{{ instructor.name }}{% if instructor.education %}, {{instructor.education}}{% endif %}</h2>
    	<h4 class="is-italic">{{ instructor.organization }}</h4>
    	<h4>{{ instructor.position }}</h4>
    	<p>{{ instructor.content | markdownify }}{% if instructor.profile_link %}<a href="{{instructor.profile_link}}" target="_blank">{{instructor.name}}'s {{instructor.profile_type}} profile</a> {% endif %}</p>
    </div>
  <br><br>
  {% endfor %}
  </section>
