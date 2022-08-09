---
permalink: /instructors
layout: page
title: Instructors
nav_order: 4
---
<!----- <h1>Instructors</h1> ---->
<!--
  {% for author in site.authors %}
  <section>
	<figure style="float:left" class="image is-128x128">
		<img class="is-rounded" src="{{ site.baseurl }}/assets/img/{{ author.picture }}">
	</figure>
    <div style="margin-left: 200px">
    	<h2>{{ author.name }}</h2>
    	<h4>{{ author.organization }}</h4>
    	<h4>{{ author.position }}</h4>
    	<p>{{ author.content | markdownify }}</p>
      {% if author.profile_link %}<br><a href="{{author.profile_link}}">{{author.short_name}}'s profile</a> {% endif %}
    </div>
  <br><br>
  {% endfor %}
  </section>
-->
  <section class="px-4">
  {% for author in site.authors %}
	<figure class="image is-128x128" style="float:left" >
		<img class="is-rounded" src="{{ site.baseurl }}/assets/img/{{ author.picture }}">
	</figure>
    <div class="px-4" style="overflow: auto">
    	<h2>{{ author.name }}</h2>
    	<h4>{{ author.organization }}</h4>
    	<h4>{{ author.position }}</h4>
    	<p>{{ author.content | markdownify }}</p>
      {% if author.profile_link %}<br><a href="{{author.profile_link}}">{{author.short_name}}'s profile</a> {% endif %}
    </div>
  <br><br>
  {% endfor %}
  </section>
