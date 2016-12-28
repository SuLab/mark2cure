SELECT  ANY_VALUE(`combined_documents`.`document_pk`) as `document_pk`,
        `combined_documents`.`pmid`,
        ANY_VALUE(`combined_documents`.`title`) as `title`,

        ANY_VALUE(`combined_documents`.`total_document_relationships`) as `total_document_relationships`,
        ANY_VALUE(`combined_documents`.`user_document_relationships`) as `user_document_relationships`,

        ANY_VALUE(`combined_documents`.`community_answered`) as `community_answered`,
        ANY_VALUE(`combined_documents`.`community_completed`) as `community_completed`,
        ANY_VALUE(`combined_documents`.`community_progress`) as `community_progress`,

        ANY_VALUE(`combined_documents`.`user_completed`) as `user_completed`,
        ANY_VALUE(`combined_documents`.`user_progress`) as `user_progress`,
        ANY_VALUE(`combined_documents`.`user_answered`) as `user_answered`,
        ANY_VALUE(`combined_documents`.`user_view_completed`) as `user_view_completed`

FROM (

  SELECT  `document_relationships`.`document_pk`,
          `document_relationships`.`pmid`,
          `document_relationships`.`title`,

          /* Note: COUNT(DISTINCT `document_relationships`.`relation_id`)
            Doesn't work as it's the total unique relationsihps that have preexisting relationship anntoations
          */
          (   SELECT COUNT(*)
            FROM `relation_relation`
            WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`
          ) as`total_document_relationships`,

          /*  This might have < 0.0 but I don't think we use it anywhere (they get filtered out with
              the user_completed boolean) and I don't honestly want to put this whole thing in an
              IF over limit statement */
          IF( (SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`) > @user_work_max,
            /* If there were more than 20 relationships, the user remaining is on those 20 */
            @user_work_max - SUM(`document_relationships`.`user_completed`),

            /* If there were 20 or less relationships, the user remaining is on them all */
            (SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`) - SUM(`document_relationships`.`user_completed`)
          ) as `user_document_relationships`,

          /* Community scope */
          SUM(`document_relationships`.`relation_unique_submissions`) as `community_answered`,
          IF(
            SUM(`document_relationships`.`relation_unique_submissions`) = (SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`)*@k_max
          , TRUE, FALSE) as `community_completed`,
          /*  This is protected from >1 b/c of the K limit per unique
              contributer limit we in the subquery */
          SUM(`document_relationships`.`relation_unique_submissions`)/((SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`)*@k_max) as `community_progress`,

          /* User specific scope */
          IF( (SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`) > @user_work_max,
            /* If there were more than 20 relationships, the user only had to do 20 */
            IF(SUM(`document_relationships`.`user_completed`) >= @user_work_max, TRUE, FALSE),

            /* If there were 20 or less relationships, the user only had to do them all */
            IF(SUM(`document_relationships`.`user_completed`) = (SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`), TRUE, FALSE)
          ) as `user_completed`,

          /*  This might have > 1.0 but I don't think we use it anywhere and
              I don't honestly want to put this whole thing in an IF over limit statement */
          IF( (SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`) > @user_work_max,
            /* If there were more than 20 relationships, the user progress is on those 20 */
            SUM(`document_relationships`.`user_completed`)/@user_work_max,

            /* If there were 20 or less relationships, the user progress is on them all */
            SUM(`document_relationships`.`user_completed`)/(SELECT COUNT(*) FROM `relation_relation` WHERE `relation_relation`.`document_id` = `document_relationships`.`document_pk`)
          ) as `user_progress`,

          SUM(`document_relationships`.`user_completed`) as `user_answered`,
          COALESCE((
            /*  Filtering out uncompleted relationship work
                that was flagged as being a completed view because
                the logic (prevous + existing) can be messy and
                we NEVER want a user to be rewarded a "Relationship
                Document Quest" bonus more than once

                I coalesce on a subquery with the `document_view`.`completed`
                b/c it's cleaner to me than joining and filtering views that
                may not exists (NULL) and uncompleted views */
              SELECT `document_view`.`completed`
              FROM `document_view`
              INNER JOIN `document_section`
                  ON `document_section`.`id` = `document_view`.`section_id`
              WHERE `document_section`.`document_id` = `document_relationships`.`document_pk`
                  AND `document_view`.`task_type` = 'ri'
                  AND `document_view`.`user_id` = @user_id
                  AND `document_view`.`completed` = True

              /*  There may be multiple Views, order first
                  by the completed ones so the LIMIT section
                  is on that */
              ORDER BY `document_view`.`completed` DESC
              LIMIT 1
          ), FALSE) as `user_view_completed`

  FROM (

    SELECT  `relationship_document`.`id` as `document_pk`,
            `relationship_document`.`document_id` as `pmid`,
            `relationship_document`.`title` as `title`,

            `relation_relation`.`id` as `relation_id`,

            /*  Counting the unique submissions on a relationship level */
            CASE
                WHEN COUNT(DISTINCT `document_view`.`user_id`) > @k_max THEN @k_max
                ELSE COUNT(DISTINCT `document_view`.`user_id`)
            END as `relation_unique_submissions`,

            /*  If the provided users of interest has provided
                an answer for the specific relationship */
            IF(
              SUM(`document_view`.`user_id` = @user_id)
            , TRUE, FALSE) as `user_completed`

    FROM (
        /*  Return back the PK, PMID and Title for all Documents
            that are in active RelationGroups */
        SELECT  `document_document`.`id` as `id`,
                `document_document`.`document_id`,
                `document_document`.`title` as `title`

        FROM `document_document`

        INNER JOIN `relation_relationgroup_documents`
            ON `relation_relationgroup_documents`.`document_id` = `document_document`.`id`

        INNER JOIN `relation_relationgroup` as `group`
            ON `group`.`id` = `relation_relationgroup_documents`.`relationgroup_id`

        WHERE `group`.`enabled` = 1

    ) as `relationship_document`

    /*  This INNER join is somewhat redundant as enforces only
        Documents with known potential relationships to show up.
        However, I'm still keeping the subquery as I think there is value
        to selectively fetch documents from a RelationGroup that have
        it's own respective controls (eg. Enabled) */
    INNER JOIN `relation_relation`
        ON `relation_relation`.`document_id` = `relationship_document`.`id`

    INNER JOIN `relation_relationannotation`
        ON `relation_relationannotation`.`relation_id` = `relation_relation`.`id`

    INNER JOIN `document_annotation`
        ON (`document_annotation`.`object_id` = `relation_relationannotation`.`id`
          AND `document_annotation`.`content_type_id` = @rel_ann_content_type_id)

    INNER JOIN `document_view`
        ON `document_view`.`id` = `document_annotation`.`view_id`

    GROUP BY `relation_id`

  ) as `document_relationships`

  GROUP BY `document_relationships`.`document_pk`

  UNION ALL

  SELECT  `document_relationships`.`id` as `document_pk`,
          `document_relationships`.`document_id` as `pmid`,
          `document_relationships`.`title` as `title`,
          (
            SELECT COUNT(*)
            FROM `relation_relation`
            WHERE `relation_relation`.`document_id` = `document_relationships`.`id`
          )  as `total_document_relationships`,
          (
            /* (TODO) Think about this field. Duplicate subquery + is it really needed? */
            SELECT COUNT(*)
            FROM `relation_relation`
            WHERE `relation_relation`.`document_id` = `document_relationships`.`id`
          )   as `user_document_relationships`,
          0 as `community_answered`,
          0 as `community_completed`,
          0.0 as `community_progress`,

          0 as `user_completed`,
        0.0 as `user_progress`,
          0 as `user_answered`,
          0 as `user_view_completed`

  FROM (
        /*  Return back the PK, PMID and Title for all Documents
            that are in active RelationGroups */
        SELECT  `document_document`.`id` as `id`,
                `document_document`.`document_id`,
                `document_document`.`title` as `title`

        FROM `document_document`

        INNER JOIN `relation_relationgroup_documents`
            ON `relation_relationgroup_documents`.`document_id` = `document_document`.`id`

        INNER JOIN `relation_relationgroup` as `group`
            ON `group`.`id` = `relation_relationgroup_documents`.`relationgroup_id`

        WHERE `group`.`enabled` = 1

  ) as `document_relationships`

) as `combined_documents`

/* (TODO) This doesn't explicitly state to take the first table from the union's results */
GROUP BY `combined_documents`.`pmid`

HAVING `community_completed` = FALSE
    AND NOT `user_completed` = TRUE
    AND NOT `user_view_completed` = TRUE
    AND NOT `total_document_relationships` = 0

ORDER BY  `community_progress` DESC,
          `user_document_relationships` DESC

LIMIT 25;
