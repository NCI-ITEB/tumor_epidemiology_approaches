---
permalink: /speakers
layout: page
title: Invited Speakers
---
<!----- <h1>Instructors</h1> ---->
<!--
  {% for speaker in site.speakers %}
  <section>
	<figure style="float:left" class="image is-128x128">
		<img class="is-rounded" src="{{ site.baseurl }}/assets/img/{{ speaker.picture }}">
	</figure>
    <div style="margin-left: 200px">
    	<h2>{{ speaker.name }}</h2>
    	<h4>{{ speaker.organization }}</h4>
    	<h4>{{ speaker.position }}</h4>
    	<p>{{ speaker.content | markdownify }}</p>
      {% if speaker.profile_link %}<br><a href="{{speaker.profile_link}}">{{speaker.short_name}}'s profile</a> {% endif %}
    </div>
  <br><br>
  {% endfor %}
  </section>
-->
  <section class="px-4">
  {% for speaker in site.speakers %}
  {% if speaker.picture %}
  <figure class="image is-128x128" style="float:left" >
		<img class="is-rounded" src="{{ site.baseurl }}/assets/speakers/{{ speaker.picture }}">
	</figure>
  {% endif %}
  <div class="px-4" style="overflow: auto">
    	<h2>{{speaker.name}}</h2>
    	<h4>{{ speaker.talk_date }} - {{ speaker.title }}</h4>
    	<p>{{ speaker.content | markdownify }}{% if speaker.profile_link %}<a href="{{speaker.profile_link}}" target="_blank">{{speaker.short_name}}'s profile</a> {% endif %}</p>
  </div>
  <br><br>
  {% endfor %}
  </section>
