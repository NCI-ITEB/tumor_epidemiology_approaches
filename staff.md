---
permalink: /instructors
layout: page
title: Instructors
---
  <section class="px-4">
  {% for author in site.authors %}
	<figure class="image is-128x128" style="float:left; clear:left" >
		<img class="is-rounded" src="{{ site.baseurl }}/assets/instructors/{{ author.picture }}">
	</figure>
    <div class="px-4" style="overflow: auto">
    	<h2>{{ author.name }}</h2>
    	<h4>{{ author.organization }}</h4>
    	<h4>{{ author.position }}</h4>
    	<p>{{ author.content | markdownify }}{% if author.profile_link %}<a href="{{author.profile_link}}" target="_blank">{{author.short_name}}'s profile</a> {% endif %}</p>
    </div>
  <br><br>
  {% endfor %}
  </section>
